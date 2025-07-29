This program is used for imputing missing covariates by the
'sequential BART' approach. The codes are revised based on the
'parallel BART' code from http://www.rob-mcculloch.org/ (though
parallel computation of BART was not used). The codes were tested in
the Linux systems (Ubuntu).

Steps in setting up the program:
1. change the working directory 'MPIDIR=' in the Makefile
2. 'cd' to the directory where Makefile was saved and run Makefile


The R function 'BartMI' is used to run the 'sequential BART for imputation'. 

R program 'example.r' illustrates how to use the function. Please
change the directory in 'source()' before using example.r.

BartMI=function(cstem,xx,yy,datatype,type,
                 numskip=199,burn=1000,
                 m=200,sigdf=3, sigquant=.90,
                 kfac=2.0)

Arguments:
1. cstem: working directory.
2. xx: covariate matrix with missing values (NAs).
3. yy: reponse (fully observed).
4. datatype: a vector indicating the type of covariates (0=continuous, 1=binary).
5. type: 0=no reponse, 1=continuous response (linear regression used for imputation) and 2=binary response (logistic regression used for imputation) 
6. numskip: number of iterations skipped
7. burn: number of iterations for burn-in
The remaining arguments are BART parameters.


    
