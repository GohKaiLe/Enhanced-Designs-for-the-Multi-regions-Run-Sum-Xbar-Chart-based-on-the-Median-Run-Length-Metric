Proc IML;

n=9;
mrl0=250;
deltamin=0.1;
deltamax=2;


Start calcR3;
rows=100;
R=J(rows,18,0);
S=J(rows,2,0);
S[1,1]=0;
S[1,2]=0;
u=1;

p2=CDF('NORMAL',3*K-delta*sqrt(n),0,1)-CDF('NORMAL',2*K-delta*sqrt(n),0,1);
p1=CDF('NORMAL',2*K-delta*sqrt(n),0,1)-CDF('NORMAL',K-delta*sqrt(n),0,1);
p0=CDF('NORMAL',K-delta*sqrt(n),0,1)-CDF('NORMAL',0-delta*sqrt(n),0,1);
p0n=CDF('NORMAL',0-delta*sqrt(n),0,1)-CDF('NORMAL',-K-delta*sqrt(n),0,1);
p1n=CDF('NORMAL',-K-delta*sqrt(n),0,1)-CDF('NORMAL',-2*K-delta*sqrt(n),0,1);
p2n=CDF('NORMAL',-2*K-delta*sqrt(n),0,1)-CDF('NORMAL',-3*K-delta*sqrt(n),0,1);

   Do a = 1 to rows;
     R[a,1]=S[a,1];   
     R[a,2]=S[a,2];
	  i=3;
	   Do j=3 to 10;
	    If (mod(j,2)=1) then do;
		                          R[a,j]=0;
							   end;
						  else do;
						          R[a,j]=R[a,2]+beta[j-i,1];
								  i=i+1;
						       end;

	   End;
      i=10;
       Do j=11 to 18;
	    If (mod(j,2)=0) then do;
		                          R[a,j]=0;
							   end;
						  else do;
						          R[a,j]=R[a,1]+alpha[j-i,1];
								  i=i+1;
						       end;

	   End;
        Do j=3 to 10;
	    If (mod(j,2)=0) then do;
		                         If (R[a,j] > beta4) then do;
								                      ww=0;
                                                      Do w = 1 to rows;
													     If (R[a,j]^=S[w,2]) then ww=ww+1;
													  End;
														 If (ww=rows) then do;
														               S[u+1,1]=0;
																	   S[u+1,2]=R[a,j];
																	   u=u+1;
																	   ww=0;
														                end;												    
													 end;
								end;
		End;
		Do j=11 to 18;
	    If (mod(j,2)=1) then do;
		                         If (R[a,j] < alpha4) then do;
								                      ww=0;
                                                      Do w = 1 to rows;
													     If (R[a,j]^=S[w,1]) then ww=ww+1;
													   End;
														 If (ww=rows) then do;
														               S[u+1,1]=R[a,j];
																	   S[u+1,2]=0;
																	   u=u+1;
																	   ww=0;
														                end;												 
													 end;
								end;
		End;
	End;

    M=J(u,rows,0);
	Do i=1 to u;
	   M[i,i]=1;
	End;

	R1=M*R;
    temp1=J(1,18,0);
    Do t=2 to u-1;
       Do v=t+1 to u;
          If(R1[t,2]>R1[v,2]) then do;		                          
		                            temp=R1[t,2];
                                    R1[t,2]=R1[v,2];
									R1[v,2]=temp;
                                    
                                    Do i=1 to 18;
                                       If(i^=2) then do;
                                          temp1[1,i]=R1[t,i];
                                          R1[t,i]=R1[v,i];
                                          R1[v,i]=temp1[1,i];
                                        end;
								     End;
								   end;
	   End;
	End;
	
    Do t=2 to u-1;
       Do v=t+1 to u;
          If(R1[t,1]>R1[v,1]) then do;		                          
		                            temp=R1[t,1];
                                    R1[t,1]=R1[v,1];
									R1[v,1]=temp;
                                    
                                    Do i=1 to 18;
                                       If(i^=1) then do;
                                          temp1[1,i]=R1[t,i];
                                          R1[t,i]=R1[v,i];
                                          R1[v,i]=temp1[1,i];
                                        end;
								     End;
								   end;
	   End;
	End;
 
	R2=J(u,9,0);
	Do b = 1 to u;
	  c=1;
	  Do a = 3 to 18 by 2;
	    d=0;
	    do row = 1 to u;
		  If (R1[b,a]=R1[row,1] & R1[b,a+1]=R1[row,2]) then do;
		                                                    R2[b,a-c]=row;
															d=1;
	                                                        end;                                                            
         end; 
         If(d^=1) then R2[b,a-c]=u+1;c=c+1; 
	   End;
	End;

	Do a=1 to u;
	   R2[a,1]=a;
	End; 

   R3=J(u,u,0);
   Pr=J(1,6,0);
   Pr[1,1]=p2n;
   Pr[1,2]=p1n;
   Pr[1,3]=p0n;
   Pr[1,4]=p0;
   Pr[1,5]=p1;
   Pr[1,6]=p2;

   Do a=1 to u;
      Do f=1 to u;
	     prob=0;
		Do b=3 to 8;
		   If R2[a,b]=f then prob=prob + Pr[1,b-2];
		End;
	       R3[a,f]=prob;
	  End;
   End;
finish calcR3;

start calcprob; 
e1=J(u,1,0);
e1[1,]=1;
one=J(u,1,1);
Prob1=e1`*(i(u)-R3**mrl0)*one;
Prob2=e1`*(i(u)-R3**(mrl0-1))*one;
finish calcprob;

start calcPercentile;
Do perc=pini to 10000; 
e1=J(u,1,0);
e1[1,]=1;
one=J(u,1,1);
Prob=e1`*(i(u)-R3**perc)*one;
If (Prob>probtemp) then do;
			percProb=perc;
			perc=10001;
			end;
End;
If(perc=10002) then perc=percProb;
finish calcPercentile;

start converg;
If abs(Prob1-0.5)<0.00001 & Prob2<=0.5 then ok=1;
else ok=0;
finish converg;

start searchK;
	Do aa=1 to 3.5 by 0.5;
		x1=aa;
		K=x1;
		run calcR3;
		run calcprob;
		f1=Prob1;

		run converg;
		If ok=1 then goto lend;

		If (f1<0.5) then do;
			aa=4;
			step=-0.1;
		end;
		If aa=3.5 then do;
			K=0;
			goto lend;
		end;
	End;

	/*bracket solution*/
    Do while (ok=0);
        x2=x1+step;
		K=x2;
		run calcR3;
		run calcprob;
		f2=Prob1;

        run converg;
	    If ok=1 then goto lend;

        ftau=(f1-0.5)*(f2-0.5);
        If ftau>0 then do;
                step=step*1.2;
                x1=x2;
                f1=f2;
                ok=0;
            end;
        else ok=1;
    end; /*end of do while loop*/

    /*swap x1 with x2 and f1 with f2*/
    If x1>x2 then do;
            tempx=x1;
            x1=x2;
            x2=tempx;

            tempf=f1;
            f1=f2;
            f2=tempf;
        end;

    /*finalize solution by halving method*/
    Do while (1);
        /*find a candidate value for x by linear interpolating x1 and x2*/
        a1=(f2-f1)/(x2-x1);
        b1=f1-a1*x1;
        x=(0.5-b1)/a1;
		K=x;
		run calcR3;
		run calcprob;
		f=Prob1;

		run converg;
	    If ok=1 then goto lend;
	
        
        /*update x1 and x2 by x*/        
        If f<0.5 then do;
                x2=x;
                f2=f;
            end;
        else do;
                x1=x;
                f1=f;
            end;

    end; /*end of do while loop*/
lend:
finish searchK;

ivalue=1;
EMRLt=1E10;

/*Abscissas for Legendre-Gauss Quadrature*/
        xij={-0.9681602395
              -0.8360311073
              -0.6133714327
              -0.3242534234
              0.0
              0.3242534234
              0.6133714327
              0.8360311073
              0.9681602395};

         /*Weight for Legendre-Gauss Quadrature*/
        weij={0.0812743884
                0.1806481607
                0.2606106964
                0.3123470770
                0.3302393550
                0.3123470770
                0.2606106964
                0.1806481607
                0.0812743884};
Do alpha1=0 to 1;
   Do alpha2=alpha1 to 3;
      Do alpha3=alpha2 to 5;
	     Do alpha4=alpha3 to 10;

beta1=-alpha1;beta2=-alpha2;beta3=-alpha3;beta4=-alpha4;
beta=J(4,1,0);
beta[1,]=beta4;
beta[2,]=beta3;
beta[3,]=beta2;
beta[4,]=beta1;

alpha=J(4,1,0);
alpha[1,]=alpha1;
alpha[2,]=alpha2;
alpha[3,]=alpha3;
alpha[4,]=alpha4;
delta=0;
run searchK;
If (K^=0) then do;
	EMRLsum=0;
       Do ij=1 to 9;
           delta=((deltamax-deltamin)/2)*xij[ij]+((deltamax+deltamin)/2); 
		   run calcR3;
		   pini=1;
		   probtemp=0.5;
		   run calcPercentile;
		   MRL1=perc;
	 EMRL=weij[ij]*(1/(deltamax-deltamin))*MRL1;
     EMRLsum=EMRLsum+EMRL;
     end;
        EMRL=((deltamax-deltamin)/2)*EMRLsum;
end;
else EMRL=1E5;

if EMRL<EMRLt then do;
	EMRLt=EMRL;
	score1=alpha1;
	score2=alpha2;
	score3=alpha3;
	score4=alpha4;
	Kf=round(K,0.0001);
	EMRLf=round(EMRL,0.01);
end;

if ivalue=1 then do;
				_EMRL=EMRL;		
                _K=K;
				_alpha1=alpha1;
                _alpha2=alpha2;
                _alpha3=alpha3;
                _alpha4=alpha4;	
			end;
			else do;	
				_EMRL=_EMRL//EMRL;
                _K=_K//K;
				_alpha1=_alpha1//alpha1;
                _alpha2=_alpha2//alpha2;
                _alpha3=_alpha3//alpha3;
                _alpha4=_alpha4//alpha4;
			end;
		ivalue=ivalue+1;
End; /*end of alpha4 loop*/
End; /*end of alpha3 loop*/
End; /*end of alpha2 loop*/
End; /*end of alpha1 loop*/

delta=0;
alpha1=score1;
alpha2=score2;
alpha3=score3;
alpha4=score4;
K=Kf;
beta1=-alpha1;beta2=-alpha2;beta3=-alpha3;beta4=-alpha4;
beta=J(4,1,0);
beta[1,]=beta4;
beta[2,]=beta3;
beta[3,]=beta2;
beta[4,]=beta1;
alpha=J(4,1,0);
alpha[1,]=alpha1;
alpha[2,]=alpha2;
alpha[3,]=alpha3;
alpha[4,]=alpha4;
run calcR3;
pini=mrl0-5;
probtemp=0.5;
run calcPercentile;
MRL0f=perc;

print score1 score2 score3 score4 Kf MRL0f EMRLf;
create EMRL var{_EMRL _K _alpha1 _alpha2 _alpha3 _alpha4};
append;

quit;

/*This part gives EMRL by a, alpha1, alpha2, alpha3 & alpha4*/
ods trace on /listing;
proc means data=EMRL maxdec=4 mean;
var _EMRL;	
class _alpha1 _alpha2 _alpha3 _alpha4 _K;	
ods output Summary=EMRL_output;
run;
ods trace off;

/*This part gives the minimum EMRL*/ 
proc means data=EMRL_output maxdec=4 min;
	var _EMRL_Mean;	
	ods output Summary=minEMRL_output;
run;
 
/*Use CALL SYMPUT to create macro variable*/
data _null_;
  set minEMRL_output;
  call SYMPUT('minEMRL',_EMRL_Mean_min);
run;

%put &minEMRL;

/*Use macro variable in subsetting WHERE clause*/
data corresponding;
  set EMRL_output(where=(round(_EMRL_Mean,.0001) = round(&minEMRL ,.0001)));
run;

/*This part prints the optimal (a, alpha1, alpha2, alpha3 & alpha4)*/ 
PROC PRINT DATA=corresponding;
var _alpha1 _alpha2 _alpha3 _alpha4 _K _EMRL_Mean; 
WHERE round(_EMRL_Mean,.01)=round(&minEMRL,.01);
RUN; 
