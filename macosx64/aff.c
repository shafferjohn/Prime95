#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#include <pthread.h>
#include <stdint.h>
#include <mach/thread_policy.h>

/* Alternatively, this can be simplified to set-only without checking core count */
int32_t pu_per_core = 0;
int32_t ncore = 0;

__attribute__((constructor))
static inline void initialize_pu_per_core() {
    int64_t cacheconfig[32];
    int32_t *cacheconfig32 = cacheconfig;
    size_t size = sizeof(cacheconfig);
    sysctlbyname("hw.cacheconfig", &cacheconfig, &size, NULL, 0);
    /* per hwloc, there are 32 and 64-bit configs... */
    if (cacheconfig[0] <= 0xFFFFFFFFUL) {
        pu_per_core = cacheconfig[2];
    } else {
        pu_per_core = cacheconfig32[2];
    }
}

__attribute__((constructor))
static inline void initialize_ncore() {
    size_t size = sizeof(ncore);
    sysctlbyname("hw.ncpu", &ncore, &size, NULL, 0);
    // assert: 0 == ncore % pu_per_core
    ncore /= pu_per_core;
}

/* Index is 1-based! 0 means no binding. Return value is positive-errno. Most commonly you get EPFNOSUPPORT (46) on ARM */
int mach_set_thread_cpubind(thread_t thread, int cpu) {
    if (cpu < 0 || cpu > ncore) {
        return EINVAL;
    }

    thread_act_t handle = pthread_mach_thread_np(thread);
    thread_affinity_policy_data_t policy = {.affinity_tag = cpu};
    return thread_policy_set(handle, THREAD_AFFINITY_POLICY, (thread_policy_t)&policy, THREAD_AFFINITY_POLICY_COUNT);
}

/* Return value is tag or negative-errno */
int mach_get_thread_cpubind(thread_t thread) {W
    thread_act_t handle = pthread_mach_thread_np(thread);
    thread_affinity_policy_data_t policy;
    mach_msg_type_number_t count = THREAD_AFFINITY_POLICY_COUNT;
    kern_return_t kr = thread_policy_get(handle, THREAD_AFFINITY_POLICY, (thread_policy_t)&policy, &count, 0);
    if (kr != KERN_SUCCESS) {
        return -kr;
    }

    return policy.affinity_tag;
}
#endif
