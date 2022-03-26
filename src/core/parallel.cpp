#include <mutex>
#include <thread>

#include "api.h"
#include "parallel.h"

CPBRT_THREAD_LOCAL int ThreadIndex;

// Local definitions.

class ParallelForLoop {
public:
    std::function<void(int64_t)> func1D;
    std::function<void(Point2i)> func2D;
    const int64_t maxIndex;
    const int chunkSize;
    const uint64_t profilerState;
    // Next iteration index to be claimed by an idle worker.
    int64_t nextIndex = 0;
    int activeWorkers = 0;
    // Linked list link. Head of the list is static workList.
    ParallelForLoop *next = nullptr;
    int nX = -1;

    ParallelForLoop(
        std::function<void(int)> func1D,
        int64_t maxIndex,
        int chunkSize,
        uint64_t profilerState
    ) : func1D(std::move(func1D)),
        maxIndex(maxIndex),
        chunkSize(chunkSize),
        profilerState(profilerState)
    {}

    ParallelForLoop(
        const std::function<void(Point2i)> &func2D,
        const Point2i &count,
        uint64_t profilerState
    ) : func2D(func2D),
        maxIndex(count.x * count.y),
        chunkSize(1),
        profilerState(profilerState) {

        nX = count.x;
    }

    bool Finished() const {
        return nextIndex >= maxIndex && activeWorkers == 0;
    }
};

static std::vector<std::thread> threads;
static bool shutdownThreads = false;

// A list of parallel loops that haven't finished. New loops are added always added
// to the head of the list; this allows for nested loops: the full enclosed loop will be
// executed before the enclosing iteration.
static ParallelForLoop *workList = nullptr;
static std::mutex workListMutex;

// Worker threads wait on this condition variable to be notified of the addition of new
// loops to the workList after having seen it empty.
static std::condition_variable workListCondition;

static void workerThreadFunc(int tIndex) {
    ThreadIndex = tIndex;

    // Creation of the unique_lock also locks the mutex.
    std::unique_lock<std::mutex> lock(workListMutex);

    while (!shutdownThreads) {
        // !workList is the condition associated with workListCondition.
        if (!workList) {
            // There aren't any loops to execute at the moment. Sleep until there are
            // more loops to run. When the condition variable is signalled by the thread
            // that added a new loop to the workList, all threads wake up and compete for
            // the input mutex (workListMutex). The  one that acquires it will resume
            // execution here.
            workListCondition.wait(lock);
        } else {
            // Get work from workList and run loop iterations.
            ParallelForLoop &loop = *workList;

            // The chunkSize determines the number of consecutive iterations that this thread
            // will execute.
            int64_t indexStart = loop.nextIndex;
            int64_t indexEnd = std::min(indexStart + loop.chunkSize, loop.maxIndex);
            
            // Update loop to reflect iterations this thread will run.
            loop.nextIndex = indexEnd;
            if (loop.nextIndex == loop.maxIndex) {
                // This thread is the one that'll execute the last chunk of the loop at the head
                // of the list, completing it. Move the head of the list to the next loop. Note
                // that the assigned value will be null if there are no more loops in the list;
                // this is the condition that tells other threads that there isn't any work to do
                // at the moment. 
                workList = loop.next;
            }
            loop.activeWorkers++;

            // Run loop indices in [indexStart, indexEnd). Unlock the mutex so that other threads
            // can proceed while this one executes the function. (Up until now, all idle worker
            // threads have been waiting on the lock, after the signalling of the condition variable
            // that woke them up.)
            lock.unlock();
            for (int64_t index = indexStart; index < indexEnd; ++index) {
                if (loop.func1D) {
                    loop.func1D(index);
                } else {
                    loop.func2D(Point2i(index % loop.nX, index / loop.nX));
                }
            }

            // Update loop to reflect completion of iterations.
            lock.lock();
            loop.activeWorkers--;
        
            if (loop.Finished()) {
                workListCondition.notify_all();
            }
        }
    }

    // TODO: Report thread statistics at worker thread exit.
}

// Definitions.

void ParallelFor(
    // Wraps a callable: a lambda expression, a function, a pointer to a function, etc.
    // The callable receives the iteration's index.
    const std::function<void(int64_t)> &func,
    // Total number of iterations.
    int64_t count,
    int chunkSize
) {
    // Run iterations immediately if not using threads or if count is small.
    // TODO: in repo, CpbrtOptions.nThreads == 1 is threads.empty() instead.
    if (CpbrtOptions.nThreads == 1 || count < chunkSize) {
        for (int64_t  i = 0; i < count; ++i) {
            func(i);
        }
        return;
    }

    // launch worker threads if needed. Only the first ParalellFor call creates the worker
    // thread pool; subsequent calls will use the existing pool; threads will persist until
    // explicitly told to shut down with ParallelCleanup().
    // TODO: in repo, thread creation takes place in ParallelInit().
    if (threads.size() == 0) {
        ThreadIndex = 0;
        // Create 1 fewer thread than there are cores: the thread executing this ParallelFor
        // call also counts.
        for (int i = 0; i < NumSystemCores() - 1; ++i) {
            threads.push_back(std::thread(workerThreadFunc, i+1));
        }
    }
    
    // Create and enqueue ParallelForLoop for this loop.
    // TODO: implement CurrentProfilerState() in stats.h.
    ParallelForLoop loop(std::move(func), count, chunkSize, 0 /*CurrentProfilerState()*/);
    workListMutex.lock();
    loop.next = workList;
    // This new loop becomes the head of the list: a nested loop that has to finish before
    // the enclosing iteration does.
    workList = &loop;
    workListMutex.unlock();
    
    // Notify worker threads of work to be done.
    // Creation of the unique_lock also locks the mutex.
    std::unique_lock<std::mutex> lock(workListMutex);
    workListCondition.notify_all();
    
    // Help out with parallel loop iterations in the current thread.
    while (!loop.Finished()) {
        // The chunkSize determines the number of consecutive iterations that this thread
        // will execute.
        int64_t indexStart = loop.nextIndex;
        int64_t indexEnd = std::min(indexStart + loop.chunkSize, loop.maxIndex);
        
        // Update loop to reflect iterations this thread will run.
        loop.nextIndex = indexEnd;
        if (loop.nextIndex == loop.maxIndex) {
            // This thread is the one that'll execute the last chunk of the loop at the head
            // of the list, completing it. Move the head of the list to the next loop. Note
            // that the assigned value will be null if there are no more loops in the list;
            // this is the condition that tells other threads that there isn't any work to do
            // at the moment. 
            workList = loop.next;
        }
        loop.activeWorkers++;

        // Run loop indices in [indexStart, indexEnd). Unlock the mutex so that other threads
        // can proceed while this one executes the function. (Up until now, all idle worker
        // threads have been waiting on the lock, after the signalling of the condition variable
        // that woke them up.)
        lock.unlock();
        for (int64_t index = indexStart; index < indexEnd; ++index) {
            if (loop.func1D) {
                loop.func1D(index);
            } else {
                loop.func2D(Point2i(index % loop.nX, index / loop.nX));
            }
        }

        // Update loop to reflect completion of iterations.
        lock.lock();
        loop.activeWorkers--;
    }
}

void ParallelFor2D(
    std::function<void(Point2i)> func,
    const Point2i &count
) {
    // Run iterations immediately if not using threads or if count is small.
    // TODO: in repo, CpbrtOptions.nThreads == 1 is threads.empty() instead.
    if (CpbrtOptions.nThreads == 1 || count.x * count.y <= 1) {
        for (int y = 0; y < count.y; ++y) {
            for (int x = 0; x < count.x; ++x) {
                func(Point2i(x, y));
            }
        }
        return;
    }

    // Launch worker threads if needed. Only the first ParalellFor call creates the worker
    // thread pool; subsequent calls will use the existing pool; threads will persist until
    // explicitly told to shut down with ParallelCleanup().
    // TODO: in repo, thread creation takes place in ParallelInit().
    if (threads.size() == 0) {
        ThreadIndex = 0;
        // Create 1 fewer thread than there are cores: the thread executing this ParallelFor
        // call also counts.
        for (int i = 0; i < NumSystemCores() - 1; ++i) {
            threads.push_back(std::thread(workerThreadFunc, i+1));
        }
    }

    // Create and enqueue ParallelForLoop for this loop.
    // TODO: implement CurrentProfilerState() in stats.h.
    ParallelForLoop loop(std::move(func), count, 0 /*CurrentProfilerState()*/);
    
    {
        // lock_guard releases the mutex when it goes out of scope.
        std::lock_guard<std::mutex> lock(workListMutex);
        loop.next = workList;
        workList = &loop;
    }

    // Notify worker threads of work to be done.
    // Creation of the unique_lock also locks the mutex.
    std::unique_lock<std::mutex> lock(workListMutex);
    workListCondition.notify_all();

    // Help out with parallel loop iterations in the current thread.
    while (!loop.Finished()) {
        int64_t indexStart = loop.nextIndex;
        // chunkSize is always 1 for 2D loops.
        int64_t indexEnd = std::min(indexStart + loop.chunkSize, loop.maxIndex);

        loop.nextIndex = indexEnd;
        if (loop.nextIndex == loop.maxIndex) {
             // This thread is the one that'll execute the last chunk of the loop at the head
            // of the list, completing it. Move the head of the list to the next loop. Note
            // that the assigned value will be null if there are no more loops in the list;
            // this is the condition that tells other threads that there isn't any work to do
            // at the moment.
            workList = loop.next;
        }
        loop.activeWorkers++;

        // Run loop indices in [indexStart, indexEnd). Unlock the mutex so that other threads
        // can proceed while this one executes the function. (Up until now, all idle worker
        // threads have been waiting on the lock, after the signalling of the condition variable
        // that woke them up.)
        lock.unlock();
        for (int64_t index = indexStart; index < indexEnd; ++index) {
            if (loop.func1D) {
                loop.func1D(index);
            } else {
                loop.func2D(Point2i(index % loop.nX, index / loop.nX));
            }
        }

        // Update loop to reflect completion of iterations.
        lock.lock();
        loop.activeWorkers--;
    }
}

void ParallelCleanup() {
    if (threads.empty()) {
        return;
    }

    {
        // lock_guard releases the mutex when it goes out of scope.
        std::lock_guard<std::mutex> lock(workListMutex);
        // The workerThreadFunc of a thread checks this value on every iteration of its main loop.
        shutdownThreads = true;
        workListCondition.notify_all();
    }

    for (std::thread &thread : threads) {
        thread.join();
    }

    threads.erase(threads.begin(), threads.end());

    shutdownThreads = false;
}

int NumSystemCores() {
    return std::max(1u, std::thread::hardware_concurrency());
}

int MaxThreadIndex() {
    return CpbrtOptions.nThreads == 0 ? NumSystemCores() : CpbrtOptions.nThreads;
}