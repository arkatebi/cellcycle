/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *       http://www.pcg-random.org
 */

/*
 * This code is derived from the full C implementation, which is in turn
 * derived from the canonical C++ PCG implementation. The C++ version
 * has many additional features and is preferable if you can use C++ in
 * your project.
 */

#include<stdio.h>
#include "pcg_basic.h"
#include<time.h>
#include<math.h>
// state for global RNGs

static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

/*
//void set_seed(int SEED, int USER_SEED)
void set_seed(int USER_SEED, int SEED)
{ 
    //if SEED is less than or equal to zero, 
    //use time*USER_SEED to initialize the random number 
    //generator, otherwise use SEED to do such initialization.
    //USER_SEED must be greater than zero for this 
    //program to function:
    //printf("%s\n","setting seed:");
    if (SEED <=0) {
        //printf("USER_SEED: %d\n", USER_SEED);
        pcg32_srandom((uint64_t) (time(NULL)*USER_SEED), 54u);
    }
    else {
        //printf("SEED: %d\n", SEED);
        pcg32_srandom((uint64_t)SEED, 54u);
    }
}
*/

void set_seed(int USER_SEED, int SEED)
{ 
    //if USER_SEED is less than or equal to zero,
    //use SEED to initialize the random number generator; 
    //otherwise, use time*USER_SEED to do such initialization.
    
    if (USER_SEED<=0) {
        pcg32_srandom((uint64_t)SEED, 54u);
    }
    else {
        //printf("USER_SEED: %d\n", USER_SEED);
        pcg32_srandom((uint64_t) (time(NULL)*USER_SEED), 54u);
    }
}

// Generate random value in (minV, maxV) following uniform distribution
double prandu(double minV, double maxV)
{
  double u;
  do {
    u = ((double)pcg32_random()/UINT32_MAX);
  } while (u==0);
  return (minV + (maxV - minV)*u); 
}

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

void pcg32_srandom(uint64_t seed, uint64_t seq)
{
    pcg32_srandom_r(&pcg32_global, seed, seq);
}

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint32_t pcg32_random()
{
    return pcg32_random_r(&pcg32_global);
}


// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = -bound % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In 
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}

uint32_t pcg32_boundedrand(uint32_t bound)
{
    return pcg32_boundedrand_r(&pcg32_global, bound);
}

/*-----------------------------------------------------------------------*/
// Generate random value in (minV, maxV) following uniform distribution
double randu(double minV, double maxV)
{
  double u;
  do {
    u = ((double)pcg32_random()/UINT32_MAX);
  } while (u==0);
  return (minV + (maxV - minV)*u); 
}



/********----- functions for stochastic modeling ----- *****************/
/*************------------START-------------****************************/
double pcg32_gasdev()	//  Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.
{
    
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;
    // if (-rng->state < 0)
    iset=0;
    if (iset == 0) {
        do {
            double u1;
            do {
                u1 = ((double)pcg32_random()/UINT32_MAX);
            } while (u1==0);
            
            v1=2.0*u1 - 1.0;
            
            do {
                u1 = ((double)pcg32_random()/UINT32_MAX);
            } while (u1==0);
            
            v2=2.0*u1 - 1.0;
            rsq=v1*v1+v2*v2;
            
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else {   iset=0;
        return gset;
    }
}

/********----- functions for stochastic modeling ----- *****************/
/*************------------END-------------*****************************/









