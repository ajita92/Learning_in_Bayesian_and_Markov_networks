
########### Input Format of training Data ###################
*.dat file contains 1000000 samples for the respective network

First line contains space separated variables of the network.
Remaining lines contain the sample values of the corresponding variable.

####### Input Format of Query File #################
The format of Query file is same as *.dat file consisting of n lines corresponding to n test points where in each line,  some variables are given values from their domain. These variables are evidence. The rest of variable are given no values denoted by “?” . The “?” Variables are Query variables.  

#####Output Format for parameters of Bayesian Network###############

You have to replace each “?” in the given BIF file with the corresponding conditional probability. This should comply with semantics of BIF file format.

####Output Format for parameters of Markov Network  ################

For Markov network, you just need to modify the probability part of bif file:

For e.g 
Bayesian Network CPT of insurance.bif is
 
probability ( SocioEcon | Age ) {
  (Adolescent) ?, ?, ?, ?;
  (Adult) ?, ?, ?, ?;
  (Senior) ?, ?, ?, ?;
}

will change into below format (Save it with .txt extension) :

probability ( SocioEcon, Age ) {
  (Prole, Adolescent) ?;
  (Prole, Adult) ?;
  (Prole, Senior) ?;
  (Middle, Adolescent) ?;
  (Middle, Adult) ?;
  (Middle, Senior) ?;	
  (UpperMiddle, Adolescent) ?;
  (UpperMiddle, Adult) ?;
  (UpplerMiddle, Senior) ?;
  (Wealthy, Adolescent) ?;
  (Wealthy, Adult) ?;
  (Wealthy, Senior) ?;
}

####### Output format of Probability Values for Query Variables #####################
Extension of the o/p of bayesian-network should be .bn.out
Extension of the o/p of markov-network should be .mn.out

Report the marginal probability of all possible assignments of a missing variable in a test case. Please use the following format to report the results:
1. All test case results should be separated by a blank line.
2. Report the prediction for every variable in a seperate line.
3. For every varaible report the marginal of all possible assignment in the following manner:
   For e.g if SocioEcon is query variable and p1-p4 are marginal probabilities for each discrete assignment (Prole, Middle, UpperMiddle, Wealthy) of the variable.
   SocioEcon Prole:p1 Middle:p2 UpperMiddle:p3 Wealthy:p4;

Sample output file:
test1-var1 test1-var1-ass1:p1 test1-var1-ass1:p2 test1-var1-ass1:p3
test1-var2 test1-var2-ass1:p1 test1-var2-ass1:p2

test2-var1 test2-var1-ass1:p1 test2-var1-ass1:p2 test2-var1-ass1:p3
test2-var2 test2-var2-ass1:p1 test2-var2-ass1:p2 test2-var2-ass1:p3
test2-var3 test2-var3-ass1:p1 test2-var3-ass1:p2 

test3-var1 test3-var1-ass1:p1 test3-var1-ass1:p2 
test3-var2 test3-var2-ass1:p1 test3-var2-ass1:p2 

####### AUC Score ##################
To Be updated



 

