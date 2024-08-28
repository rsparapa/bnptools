
proc import datafile='hnscc.csv' out=hnscc;
guessingrows=max;
run;

proc freq;
    tables p16n*p16m*pniy / list;
run;

data hnscc;
    set hnscc;
    if p16n=p16m=0 then pniy=0;
run;

proc freq;
    tables p16n*p16m*pniy / list;
run;

ods graphics on / reset imagename="diradj" imagefmt=pdf;

proc phreg;
    class study;
    model DFS_Time*DFS(0) = age White Black Asian Male Smoker Alcoholic 
        PriorChemo Larynx OralCavity Oropharynx /*P16N P16M*/ ECSY ECSN ECSM 
        Positve_MY Positve_MN Positve_MM Close_MY Close_MN Close_MM 
/*PNIY*/ PNIN PNIM LVIY LVIN LVIM LN2Y LN2N LN2M RSY RSN RSM Study 
        diag_year ;
    assess ph / seed=10;
    strata P16N P16M PNIY; 
    baseline out=surv lower=lower survival=surv upper=upper / 
        diradj group=study;
run;

/*
%_pdfjam;

%sysexec mv diradj.pdf assess.pdf;
*/

data surv;
    set surv;
    if lower=. then lower=1;
    if upper=. then upper=1;
run;

proc print;
    var DFS_Time study P16N P16M PNIY lower surv upper;
    format _numeric_ 4.2;
run;

proc export replace data=surv outfile='diradj.csv';
run;

footnote;

proc sgpanel data=surv;
    panelby P16N P16M PNIY;
    series x=dfs_time y=surv / group=study;
run;

/*
%_pdfjam;

%sysexec mv diradj21.pdf diradj.pdf;
*/
