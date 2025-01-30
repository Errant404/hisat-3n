/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef WITH_TBB
# include <tbb/mutex.h>
# include <tbb/spin_mutex.h>
# include <tbb/queuing_mutex.h>
#else
# include <mutex>
#endif

#ifdef NO_SPINLOCK
# ifdef WITH_TBB
#   ifdef WITH_QUEUELOCK
#  	define MUTEX_T tbb::queuing_mutex
#   else
#       define MUTEX_T tbb::mutex
#   endif
# else
#   define MUTEX_T std::mutex
# endif
#else
# ifdef WITH_TBB
#   define MUTEX_T tbb::spin_mutex
# else
#   define MUTEX_T std::mutex
# endif
#endif /* NO_SPINLOCK */


/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:

	ThreadSafe() : ptr_mutex(NULL) { }
	
	ThreadSafe(MUTEX_T* ptr_mutex, bool locked = true) : ptr_mutex(NULL) {
		if(locked) {
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
			//have to use the heap as we can't copy
			//the scoped lock
			this->ptr_mutex = new MUTEX_T::scoped_lock(*ptr_mutex);
#else
			this->ptr_mutex = ptr_mutex;
			ptr_mutex->lock();
#endif
		}
	}

    ~ThreadSafe() {
        unlock();
    }

    void unlock() {
        if (ptr_mutex != NULL) {
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
            delete ptr_mutex;
#else
            ptr_mutex->unlock();
#endif
            ptr_mutex = NULL;
        }
    }

private:
#if WITH_TBB && NO_SPINLOCK && WITH_QUEUELOCK
	MUTEX_T::scoped_lock* ptr_mutex;
#else
	MUTEX_T *ptr_mutex;
#endif
};

#endif
