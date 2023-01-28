#include <iostream>
#include <fstream>
#include <sstream>

#include <stdio.h>
#include <vector>
#include <tuple>
#include <utility>
#include <bitset>
#include <math.h>
#include <string>
#include <chrono>

typedef unsigned long int uint32;
typedef unsigned long long int uint64;
typedef __uint128_t uint128;

std::pair<uint32,uint64> factorTwos(uint64 number){
    if(number == 1){
        return std::make_pair(0,1);
    }
    uint64 x = number;
    uint32 exp = 0;
    uint64 mult = 1;

    while((x & 1) == 0){
        x = x >> 1;
        exp++;
    }
    mult = x;

    return std::make_pair(exp,mult);
}

/*
    Use the Square-Multiply Algorithm to calculate (base^y) mod modulus
    The binary version only works on 32-bit numbers, as 64-bit numbers will overflow.
        this runs in O(log_2 y) time, since we go through the loop for each digit of y in binary.

*/
uint64 expmod(uint64 base, uint64 y, uint64 modulus){
    if(base == 0){
        return 0;
    }
    if(base == 1){
        return 1;
    }

    uint64 exp = y;

    uint128 res = 1;

    uint128 factor = base % modulus;

    while(y > 0){
        //if y is odd, multiply by factor
        if(y & 1)
            res = (res*factor)%modulus;
        //std::cout << "y: " << y << " res: " << res << std::endl;
        //now y must be even; divide by two
        y = y >> 1; //y = y/2
        factor = (factor*factor) % modulus;
        //std::cout << "factor: " << factor << std::endl;
    }

    return res;
}

/*
    return true if x = 0 (mod 2^t - 1)

    x = (x mod M) + floor(x/M) (mod M-1)
    this is just the case where M = 2^t

    2^t     = 1   0   ... 0 0
    2^t - 1 = 0   1   ... 1 1
              t (t-1) ... 1 0

    example: t = 3
    2^t = 8, N = 7
    1 << 3     = 1000
    1 << 3 - 1 = 0111
*/
uint64 mod_mersenne(uint64 x, uint64 t){
    uint64 modulus = ((x & ((1 << t) - 1)) + (x >> t)) % ((1 << t) - 1);
    return modulus;
}

bool any_mersenne(uint64 x){
    //iterate through every mersenne prime less than x, and check if it divides x.
    //this returning true means the number cannot be prime.

    /*M_t     = 0011...11
      M_(t+1) = 0111...11
      M_(t+1) = (M_t << 1) + 1

    */
    //M_1 = 2 - 1 = 1, which isn't prime
    // start at 2^2 - 1 = 3
    uint64 mask = 3;
    uint64 t = 2;
    uint64 remain;
    //shift = x >> t, which is the same as x / 2^t
    uint64 shift = (x >> 2);
    while(shift){
        remain = ((x & mask) + shift) % mask;
        if(!remain) return true;
        shift = shift >> 1;
        mask = (mask << 1) + 1;
    }
    return false;
}

/*
given p-1 factored as a power of two times an odd number:
p-1 = 2^(expo) * mult
*/
bool mr_test(uint32 expo, uint64 s_mult, uint64 n_mult, int num_witnesses, uint64* witnesses, bool verbose=false){
    // mult == 0 -> p-1 == 2^(expo)*0 == 0 -> p == 1
    if((s_mult == 0) && (n_mult == 0))
        return false;
    
    // expo == 0 and mult == 1 -> p-1 = 1 -> p == 2
    if((expo == 0) && (s_mult == 1) && (n_mult == 1))
        return true;
    // expo == 0, mult is odd -> p-1 is odd -> p is even and not two
    if(expo == 0)
        return false; 

    uint64 p =  ((s_mult*n_mult) << expo) + 1;
    // p-1 = 2^(expo) * mult
    // p   = 2^(expo) * mult + 1

    uint64 factor;
    uint64 s_factor;

    bool witness_passed = true;
    //uint64 mult = (uint64)

    //for each witness:
    for(int i = 0; i < num_witnesses && witness_passed; i++){
            if(verbose) std::cout << "checking witness: " << witnesses[i] << std::endl;
            
            witness_passed = false;
            s_factor = expmod(witnesses[i],s_mult,p);
            

            factor = expmod(s_factor,n_mult,p);
            if(verbose) std::cout << "factor = " << factor << std::endl;

            //first test succeeds, we can move on.
            if(factor == 1 || factor == p-1){
                if(verbose) std::cout << "witness = " << witnesses[i] << " certifies." << std::endl;
                witness_passed = true;
                continue;
            }

            for(uint64 i = 1; i < expo && !witness_passed; i++){
                    //square it repeatedly
                    if(verbose) std::cout << "entered for loop for exponents. i = " << i << std::endl;
                    factor = (uint64)(((uint128)factor*(uint128)factor) % p);
                    if(verbose) std::cout << "new factor = " << factor << std::endl;

                    //check a^(2^(i+1))xz = -1 (mod p)
                    // = (a^(2^i * xz)) ^ 2
                    if(factor == p-1){
                        //this witness certifies it as prime. Go to the next witness.
                        if(verbose) std::cout << "witness = " << witnesses[i] << " certifies with expo 2^" << i << std::endl;
                        witness_passed = true;
                        break;
                    }
                    if(factor == 1){
                        //this witness ensures the number is composite.
                        if(verbose) std::cout << "factor == 1 for witness = " << witnesses[i] << " at exponent: 2^" << i << std::endl; 
                        return false;
                    }

            }
    }
    if(verbose) std::cout << "reached end. witness_passed = " << witness_passed << std::endl;
    return witness_passed;


}

bool is_prime(uint64 p, int num_witnesses, uint64* witnesses){
    std::pair<uint32,uint64> p_twos = factorTwos(p-1);
    //std::cout << "exponent: " << p_twos.first << " mult: " << p_twos.second << std::endl;
    return mr_test(p_twos.first, p_twos.second, 1, num_witnesses, witnesses);
}

/*
Go through s+1, 2s+1, ... sn+1 until we find a prime p. Once we do, return it.
*/
uint64 test_s(uint64 s_val, int num_witnesses, uint64* witnesses,bool verbose=true){
    uint64 p = 1;
    
    //factor s = 2^w * x
    std::pair<uint32,uint64> s_twos = factorTwos(s_val);
    uint32 s_exp = s_twos.first;
    uint64 s_mul = s_twos.second;

    if(verbose)
        std::cout << "S: 2^" << s_exp << " * " << s_mul << std::endl;

    if(s_exp > 63){
        std::cout << "exponent for S: " << s_exp << " is far too big to be stored." << std::endl;
    }

    uint32 n_exp = 0;
    uint64 n_mul = 1;

    std::pair<uint32,uint64> n_twos;

    uint64 n = 0;
    while(true){
        n++;
        p += s_val;
        
        //std::cout << "p: " << p << std::endl;
        //factor n = 2^y * z
        n_twos = factorTwos(n);
        n_exp = n_twos.first;
        n_mul = n_twos.second;

        //if p is even, we can immediately move on.
        //p is even -> 2^t * d + 1 is even -> 2^t * d is odd -> t = 0
        if(s_exp == 0 && n_exp == 0)
            continue;

        if(verbose)
            std::cout << "testing n = " << n << std::endl;

        if(mr_test(s_exp+n_exp, s_mul, n_mul, num_witnesses, witnesses, verbose))
            return p;
       
    } 
}

std::pair<uint64,uint64> last_spot(std::string filename ="C:\\Dirichlet_C\\exhaustive_results.csv"){

    std::ifstream myFile(filename);
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    uint64 lastS = 0;
    uint64 lastP = 0;
    std::string line;
    //std::getline(myFile,line);

    while(std::getline(myFile, line,'\n'))
    {
        std::stringstream myLine(line);
        std::string last_line;
        std::getline(myLine,last_line,',');
        lastS = std::stoull(last_line);
        std::getline(myLine,last_line,',');
        lastP = std::stoull(last_line);
    }
    myFile.close();

    return std::make_pair(lastS,lastP);
}

bool add_pair(uint64 s_val, uint64 p_val, std::string filename ="C:\\Dirichlet_C\\exhaustive_results.csv"){
    std::ofstream outfile;
    //std::cout << "added s = " << s_val << std::endl;
    outfile.open(filename, std::ios_base::app);
    if(outfile.is_open()){
        outfile << '\n' << std::to_string(s_val) << "," << std::to_string(p_val);
        outfile.close();
        return true;
    } else {
        return false;
    }
    return true;
}

void exhaustive_iter(uint64 steps, uint64* witnesses, int num_witnesses, std::string filename="C:\\Dirichlet_C\\exhaustive_results.csv", bool percentage=true){
    auto lastPlace = last_spot(filename);
    uint64 currS = lastPlace.first;
    uint64 currP = lastPlace.second;

    uint64 maxP = currP;

    int progress = -10;
    uint64 ten_per = steps/10;
        

    for(uint64 i = 0; i < steps; i++){
        if(i%ten_per == 0 && percentage){
            progress += 10;
            std::cout << progress << "%" << std::endl;

        }
        currS++;
        currP = test_s(currS,num_witnesses,witnesses,false);
        std::cout << "curr P: " << std::endl;
        if(currP > maxP){
            add_pair(currS, currP);
            std::cout << "added s = " << currS << std::endl;
            maxP = currP;
        }
    }
    std::cout << "final S: " << currS << std::endl;

}

void exhaustive_iter_start(uint64 steps, uint64* witnesses, int num_witnesses, uint64 startplace, std::string filename="C:\\Dirichlet_C\\exhaustive_results.csv", bool percentage=true){
    auto lastPlace = last_spot(filename);
    uint64 currS = startplace;
    uint64 currP = lastPlace.second;

    uint64 maxP = currP;

    int progress = -10;
    uint64 ten_per = steps/10;
        

    for(uint64 i = 0; i < steps; i++){
        if(i%ten_per == 0 && percentage){
            progress += 10;
            std::cout << progress << "%" << std::endl;

        }
        currS++;
        currP = test_s(currS,num_witnesses,witnesses,false);
        if(currP > maxP){
            add_pair(currS, currP);
            std::cout << "added s = " << currS << std::endl;
            maxP = currP;
        }
    }
    std::cout << "final S: " << currS << std::endl;

}


int main(){
    uint64 num_witnesses = 7;
    uint64 witnesses[] {2, 3, 5, 7, 11, 13, 17};

    auto place = last_spot();
    //std::cout << "last s: " << place.first << ", " << place.second << std::endl;

    /*
    std::cout << "2 power: " << expmod(2,5*31,311) << std::endl;
    std::cout << "3 power: " << expmod(3,5*31,311) << std::endl;
    std::cout << "5 power: " << expmod(5,5*31,311) << std::endl;
    
    //2^1386 = 1 (mod 1387)

    std::cout << "should get 1:" << std::endl;
    std::cout << expmod(2,1386,1387) << std::endl;
    std::cout << "should get 512:" << std::endl;
    std::cout << expmod(2,693,1387) << std::endl;
    std::cout << "should get 742:" << std::endl;
    std::cout << expmod(2,1762,1763) << std::endl;
    */
    //s =    14,317,547
    //p = 2,834,874,307
    //add_pair(255,255);
    /*
    uint64 s_to_test = 14317547;
    uint64 exp_result = 2834874307;
    uint64 n_val = (uint64)((exp_result-1)/s_to_test);
    
    uint64 result = test_s(14317547, num_witnesses, witnesses, false);
    
    std::cout << "Gotten result: " << result << std::endl;
    std::cout << "actual result: " << exp_result << std::endl;
    */

    //time test:

    //     3,894,698,099
    // 1,589,036,824,393
    
    uint64 s_to_test = 3894698099;
    uint64 exp_result = 1589036824393;

    bool expon = false;
    bool mers  = false;
    bool primality = false;
    bool timing_s = false;
    bool timing_iter = true;

    if(expon){
        //198,629,603,049
        uint64 power = 198629603049;
        uint64 mod_result = expmod(2,power,exp_result);
        std::cout << "should get 694,299,480,401" << std::endl;
        std::cout << "actually got " << mod_result << std::endl;
    }

    if(mers){
        uint64 numb = 96;
        std::cout << "testing x = " << numb << std::endl;
        std::cout << "mod 1:  " << mod_mersenne(numb, 1) << std::endl;
        std::cout << "mod 3:  " << mod_mersenne(numb, 2) << std::endl;
        std::cout << "mod 7:  " << mod_mersenne(numb, 3) << std::endl;
        std::cout << "mod 15: " << mod_mersenne(numb, 4) << std::endl;

        std::cout << "func says: " << any_mersenne(numb) << std::endl;
    }

    if(primality){
        std::cout << "check primality of: " << exp_result << std::endl;

        bool exp_is_prime = is_prime(exp_result, num_witnesses, witnesses);

        std::cout << "is it prime? " << exp_is_prime << std::endl;
        std::cout << "it should be." << std::endl;

        /*correct results:
          2^mult                           =   694,299,480,401
          2^(2*mult) = 2^(397,259,206,098) = 1,589,036,824,392
          2^(4*mult) = 2^(794,518,412,196) = 1
          
          factor   =   694,299,480,401
          factor^2 = 1,454,784,752,956
          factor^4 =   990,295,135,757
        */
    }
    
    if(timing_s){
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        uint64 result = test_s(s_to_test, num_witnesses, witnesses, false);

        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        std::chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Gotten result: " << result << std::endl;
        std::cout << "actual result: " << exp_result << std::endl;
        std::cout << diff.count() << "ms" << std::endl;
    }

    if(timing_iter){

        uint64 steps = 1000000000;
        uint64 lastended = 43639728511;

        std::cout << "starting search..." << std::endl;
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        exhaustive_iter_start(steps, witnesses, num_witnesses, lastended);

        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        std::chrono::seconds diff = std::chrono::duration_cast<std::chrono::seconds>(end - start);

        std::cout << "took " << diff.count() << "s" << std::endl;
        std::cout << "for " << steps << " iterations." << std::endl;
    }
   
}

