# HMM_algorithms
Implementing algorithms for hidden Markov models for KTH Artificial Intelligence DD2380

## HMM1
Given the current state probability distribution, return the probabilities for the 
different emissions after the next transition.

> `python hmm1.py < samples/sample_00.in`
Output: sample_00.ans

## HMM2 - forward algorithm
Calculate the probability of observing a certain emission sequence given a hidden Markov model.
> `python hmm1.py < samples/hmm2.in`
Output: hmm2.ans

## HMM3 - Viterbi algorithm
Estimate the most likely sequence of states.
> `python hmm1.py < samples/hmm3.in`
Output: hmm3.ans

## HMM4 - Baum-Welch algorithm
Estimate parameters of a hidden Markov model given an emission sequence and initial guesses for lambda.
> `g++ -std=c++11 hmm4.cpp -o hmm4`
> `./hmm4 < samples/hmm4_01.in`