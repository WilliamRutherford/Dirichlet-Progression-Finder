# Dirichlet Progression Finder
 Finds the smallest primes of the form sn+1

## Results

The results are stored in ascending order, strictly increasing. If there is a gap between these numbers, it means the result was smaller than a previous result. This is sufficient to find an upper bound, and avoids needing to store 2^32 pairs of 32-bit numbers.

Given a value s, we return a value p such that sn+1 = p is the smallest prime of this form.

This data can be used to compute a potential upper bound, and the upper bound of the data is less than any data previously published.

## Use

The function "mr_test" performs a [deterministic Miller-Rabin primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test) on a number. The larger the number, the more witness numbers required. 

"mod_mersenne" quickly checks if a mersenne prime (2^t-1) is a factor of the number.

"expmod" takes two 64-bit numbers (base and exponent) and a modulus, and returns the exponentiation. This uses the highly efficient Square-Multiply algorithm, which runs in O(log_2 y) time, where log_2 y is the number of digits of the exponent.  
