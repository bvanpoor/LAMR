//******************************************************
//	Programmer: Brett van Poorten
//	Project Name: Length-Based Population Dynamics Model
//	Date: 21 September 2011
//	Version: 1000
//	Comments: This model estimates numbers at age and uses mean growth rates from LWFitter 
//		to estimate length distributions in each sampling event
//	
//******************************************************/

// CHANGE LENGTH DISTRIBUTIONS SO LENGTHS = muLa * ( 1 + groups * cvl * 3 ) TO REPRESENT 3 STANDARD DEVIATIONS FROM THE MEAN

TOP_OF_MAIN_SECTION
  arrmblsize = 500000000;  										
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(120000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(15000000000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(2800);

DATA_SECTION
	int debug;
	!! debug=0;						//if debug=1, spits out all kinds of info
	int detail;
	!! detail=0;						//if detail=1, spits out information to help pinpoint errors
	int est;
	!! est=1;							//if est=1, estimate real parameters; otherwise simulate fake data
	int seed;	
	int stock;
	!! stock=-3;
	int outLen;
	!! outLen=0;
	LOC_CALCS
		if(est==0)	{system("Recruitment_ModelSim3.exe");	cout<<"data generated"<<endl;}					//call the data file generator before estimating anything
		ifstream ifs("seed.txt");
		ifs>>seed;
	END_CALCS

	init_int nlakes;																	//number of lakes being modelled
	init_int nbins;																	//number of length bins
	init_int delt_l;																	//width of length bins
	init_number min_l;
	number dw;																				//half of bin width
	!! dw=delt_l/2.;
	init_int ngroups;
	vector groups(1,ngroups);
	init_int ngears;
	init_number DayT;																//time interval
	init_int MaxAge;																//maximum age modelled - more than max observed
	init_number alw;																	//intercept of length-weight relationship
	init_number plw;																	//power term of length-weight relationship
// BIOENERGETICS (GROWTH) FUNCTION PARAMETERS
	init_number TySpawn;															//time of year that spawning occurs (in proportion of year)
	init_number Pgonad;															//proportion of mass lost to gonads
	init_number Wma;																//weight at maturity
	init_number Tinact;															//
	init_number Ktemp;																//
	init_number rbase;																//
	init_number rs;			
	init_number WinfFactor;														//
	init_number WsPower;															//
	init_vector Tmax(1,nlakes);													//Tmax of temperature sine function for each lake
	init_vector  Tmin(1,nlakes);													//Tmin of temperature sine function for each lake
	init_vector Twinter(1,nlakes);												//Twinter of temperature sine function for each lake
	init_vector Tyoff(1,nlakes);													//Tyoff of temperature sine function for each lake
	init_number LengthHatch;													//length at hatch
	init_number TyHatch;															//time of year when hatching occurs
	init_ivector tnyears(1,nlakes);												//number of years where tagged fish exist and where growth was estimated from growth model
	ivector glen(1,nlakes);
	!! glen=tnyears-1;
	init_imatrix gconind(1,nlakes,1,glen);										//index for Gcon's
	init_ivector datayears(1,nlakes);											//number of lake-years for which length/age data exists
	ivector gconlen(1,nlakes);
	!! gconlen=datayears-1;
	init_ivector tstart(1,nlakes);
// RECRUITMENT FUNCTION PARAMETERS
	init_int irec;
	init_ivector nyears(1,nlakes);												//number of years each lake was studied
	int tyears;																			//maximum number of years across lakes (not needed for 1-lake evaluations)
	init_ivector startyr(1,nlakes);
	init_number lr;																	//reference length for Lorenzen (2000) mortality
	init_int ObsAge;																	//number of ages seen in the catch
	ivector thetayears(1,nlakes);
	!! thetayears=nyears+ObsAge-startyr-irec+1;
	init_ivector intheta(1,nlakes);					//early recruitment isn't estimated so intheta is the number of recruits actually estimated
	int mintheta;
	!! mintheta=sum(intheta);
	init_vector A(1,nlakes);														//area of each lake
	init_int c;
	init_int ndays;																	//number of days per year being modelled (365/10)
	int yrs1;
	int yrs2;
	int yrs3;
	int yrs4;
	int yrs5;
	int yrs6;
	int yrs7;
	int yrs8;
	int studlen;
	init_int dend;
	!!if(debug==1)	cout<<"data reading done"<<endl;
	!! if(dend!=999)	{cout<<"input error! dend = "<<dend<<endl;	exit(1);}
	
	!! ad_comm::change_datafile_name("Catches.dat");
	init_ivector ncatch(1,nlakes);
	init_3darray CatchDat(1,nlakes,1,ncatch,1,12);
	init_ivector neffort(1,nlakes);
	init_3darray EffortDat(1,nlakes,1,neffort,1,7);
	int meffort;
	!! meffort=max(neffort);
	init_int dend4;
	!! if(debug==1)	cout<<"catch reading done"<<endl;
	!! if(dend4!=999)	{cout<<"catch input error! dend= "<<dend4<<endl;	exit(1);}

	int Amax;
	int ObsAmax;
	int mage;
	!! mage=ObsAge+1;
	vector bins(1,nbins);
	vector days(1,ndays);
	vector tdays(1,nlakes);
	!! tdays=(nyears)*ndays;
	vector Tbar(1,nlakes);

	matrix TyCatch(1,nlakes,1,ncatch);											//date each fish was captured as a proportion of the year
	matrix CID(1,nlakes,1,ncatch);												//ID of each fish (0=not in recap database; #=reference number in recap database)
	matrix CLen(1,nlakes,1,ncatch);												//length of each fish at capture
	matrix CMort(1,nlakes,1,ncatch);												//captured fish that died during handling
	matrix CMark(1,nlakes,1,ncatch);											//captured fish that were marked before release
	matrix CYear(1,nlakes,1,ncatch);												//year fish was captured
	matrix CRecap(1,nlakes,1,ncatch);											//captured fish that were marked earlier that year
	3darray Ct(1,nlakes,1,ngears,1,ncatch);									//binary: either caught by gear (1) or not
	matrix CAge(1,nlakes,1,ncatch);												//assessed age of individuals upon capture
	matrix EDate(1,nlakes,1,neffort);											//date of each capture event (in proportion of year)
	matrix EYear(1,nlakes,1,neffort);											//year of each capture event relative to first year
	matrix Eflag(1,nlakes,1,neffort);											//if flag==1, the next effort date was only one (or a couple) of days later, but still processed as if the same time-step (i.e. same size of fish)
	3darray Et(1,nlakes,1,ngears,1,neffort);									//number of units of effort from each of the gear-types; right now, they are in the order of 1) fyke net; 2) gill net; 3) mini-fyke net

	4darray Catch(1,nlakes,1,meffort,1,ngears,1,6);

	5darray CUnmarkFreq(1,nlakes,1,ngears,1,meffort,irec,mage,1,nbins);			//number of recaps that died in each length bin at each sampling occasion
	5darray CMarkFreq(1,nlakes,1,ngears,1,meffort,irec,mage,1,nbins);			//number of unmarks that died in each length bin at each sampling occasion

	LOC_CALCS
		if(debug==1)	cout<<"start 1st calculations"<<endl;
		for(int gt=1;gt<=ngroups;gt++)
			groups(gt)=gt-(ngroups+1)/2.;
		groups/=max(groups);
		tyears=max(nyears);
		studlen=tyears*ndays*nlakes;
		yrs1=intheta(1);
		yrs2=intheta(2);
		yrs3=intheta(3);
		yrs4=intheta(4);
		yrs5=intheta(5);
		yrs6=intheta(6);
		yrs7=intheta(7);
		yrs8=intheta(8);
		ObsAmax=ObsAge*ndays+1.;
		bins.fill_seqadd(min_l+dw,delt_l);	
		days.fill_seqadd(DayT,DayT);
		Amax=int(MaxAge*ndays)+1;
		Catch.initialize();
		CUnmarkFreq.initialize();
		CMarkFreq.initialize();
		for(int lake=1;lake<=nlakes;lake++)
		{
			TyCatch(lake)=column(CatchDat(lake),1);
			CID(lake)=column(CatchDat(lake),2);
			CLen(lake)=column(CatchDat(lake),3);
			CMort(lake)=column(CatchDat(lake),4);
			CMark(lake)=column(CatchDat(lake),5);
			CYear(lake)=column(CatchDat(lake),6);
			CRecap(lake)=column(CatchDat(lake),7);
			CAge(lake)=column(CatchDat(lake),9+ngears);
			EDate(lake)=column(EffortDat(lake),1);
			EYear(lake)=column(EffortDat(lake),2);
			Eflag(lake)=column(EffortDat(lake),ngears+3);
			for(int gear=1;gear<=ngears;gear++)
			{
				Ct(lake,gear)=column(CatchDat(lake),7+gear);
				Et(lake,gear)=column(EffortDat(lake),2+gear)/A(lake);
			}
			if(debug==1)	cout<<"lake"<<endl;
			for(int i=1;i<=ncatch(lake);i++)
			{
				int j=0;
				for(int k=1;k<=neffort(lake);k++)
				{	if(int(EDate(lake,k)*1000)==int(TyCatch(lake,i)*1000) && EYear(lake,k)==CYear(lake,i))  j=k;}
				int bin=(CLen(lake,i)+dw-min_l)/delt_l;
				if(bin>1)
				{
					for(int gear=1;gear<=ngears;gear++)
					{
						if(Ct(lake,gear,i)>0.)
						{if(detail==1)	cout<<"lake "<<lake<<"\t"<<"i "<<i<<"\t"<<j<<"\t"<<"gear\t"<<gear<<"\t"<<int(EDate(lake,j)*1000)<<" "<<int(TyCatch(lake,i)*1000)<<"\t"<<EYear(lake,j)<<"\t"<<CYear(lake,i)<<endl;
							Catch(lake,j,gear,1)+=CMark(lake,i)*(1.-CMort(lake,i));								//fish marked and returned to the lake (alive)
							Catch(lake,j,gear,2)+=CMark(lake,i)*CMort(lake,i);									//fish marked but died
							Catch(lake,j,gear,3)+=(1.-(CMark(lake,i)+CRecap(lake,i)))*CMort(lake,i);		//unmarked fish that died
							Catch(lake,j,gear,4)+=CRecap(lake,i)*CMort(lake,i);									//recaped fish that died
							Catch(lake,j,gear,5)+=CRecap(lake,i)*(1.-CMort(lake,i));								//recaped fish that lived
							Catch(lake,j,gear,6)+=(1.-(CMark(lake,i)+CRecap(lake,i)))*(1.-CMort(lake,i));	//unmarked fish that lived but weren't marked
							if(CAge(lake,i)<irec)	{
								CUnmarkFreq(lake,gear,j,mage,bin)+=(1-CRecap(lake,i));//*C1stDay(lake,i);
								CMarkFreq(lake,gear,j,mage,bin)+=CRecap(lake,i);//*C1stDay(lake,i);
							}	else	{
								int age=CAge(lake,i)+1;
								CUnmarkFreq(lake,gear,j,age,bin)+=(1-CRecap(lake,i));//*C1stDay(lake,i);
								CMarkFreq(lake,gear,j,age,bin)+=CRecap(lake,i);//*C1stDay(lake,i);
							}
						}
					}
				}
			}
			//cout<<"catch allocation done"<<endl;
		}
		if(debug==1)	cout<<"first calcs done"<<endl;
	END_CALCS

	vector Ti(1,Amax);
	matrix Temp(1,nlakes,1,365);
	number tiny;
	int endyr;
	int recage;
	int myclass;
	LOC_CALCS
		Ti.fill_seqadd(TyHatch,DayT);
		for(int lake=1;lake<=nlakes;lake++)
		{
			for(int i=1;i<=365;i++)
			{
				Temp(lake,i)=Tmin(lake)+(Tmax(lake)-Tmin(lake))*sin((i/365.+Tyoff(lake))*2.*3.14159);
				if(Twinter(lake)>Temp(lake,i))	{Temp(lake,i)=Twinter(lake);}
			}
			Tbar(lake)=mean(Temp(lake));
		}
		tiny=1e-100;
		endyr=(1.-TyHatch)*ndays;
		recage=(irec-TyHatch)*ndays;
		myclass=max(thetayears);
	END_CALCS
	!!if(debug==1)	cout<<"data section done"<<endl;
	
	int ngro_par;
	!! ngro_par=sum(datayears)+3;
	
	!! ad_comm::change_datafile_name("growth means.txt");
	init_matrix gro_mean(1,1000,1,ngro_par);
	
	!! ad_comm::change_datafile_name("grow sds.txt");
	init_matrix gro_sd(1,1000,1,ngro_par);
	
PARAMETER_SECTION
	init_vector lM(1,nlakes,4);														//instantaneous mortality
	init_matrix lq(1,nlakes,1,ngears,2);
	init_vector lqalpha(1,ngears,3);												//alpha of selectivity
	init_vector lqbeta(1,ngears,3);													//beta of selectivity
	init_bounded_vector tqgamma(1,2,-100.,15.,3);					//gamma of selectivity
	init_bounded_number lCVl(-10.,5.,3);													//cv of length-at-age
	init_matrix ltheta(1,nlakes,1,intheta,-1);
	init_vector mutheta(1,nlakes,1);
	init_vector lsigtheta(1,nlakes,1);
	init_bounded_dev_vector Ctheta(1,yrs1,-10.,10.,1);
	init_bounded_dev_vector Caththeta(1,yrs2,-10.,10.,1);
	init_bounded_dev_vector Dtheta(1,yrs3,-10.,10.,1);
	init_bounded_dev_vector Megtheta(1,yrs4,-10.,10.,1);
	init_bounded_dev_vector Mtheta(1,yrs5,-10.,10.,1);
	init_bounded_dev_vector Ntheta(1,yrs6,-10.,10.,1);
	init_bounded_dev_vector Ptheta(1,yrs7,-10.,10.,1);
	init_bounded_dev_vector Wtheta(1,yrs8,-10.,10.,1);
	init_number Mcon(-1);
	init_number Gpow(-1);
	init_number Mpow(-1);
	init_number Qten(-1);
	init_number QtenCon(-1);
	init_vector Gconinit(1,nlakes,-1);
	init_matrix Gconmults(1,nlakes,1,gconlen,-1);
	init_number dend4(-1);
	sdreport_vector logD(1,mintheta);												//log-densities to use for plotting
	sdreport_vector preddens(1,studlen);								//predicted density of fish across lakes throughout the study
	
	objective_function_value f;

	vector pM(1,nlakes);
	matrix pq(1,nlakes,1,ngears);
	matrix palpha(1,nlakes,1,ngears);
	matrix pbeta(1,nlakes,1,ngears);	
	matrix pgamma(1,nlakes,1,ngears);
	vector pCVl(1,nlakes);
	vector sigtheta(1,nlakes);
	matrix rtheta(1,nlakes,1,intheta);
	matrix ptheta(1,nlakes,1,thetayears);
	matrix Gcon(1,nlakes,1,tnyears);
	3darray Lage(1,nlakes,1,Amax,1,nyears);
//	3darray count(1,meffort,1,2,1,nbins);

	4darray PltU(1,nlakes,1,meffort,1,ngears,1,nbins);			//number of unmarks captured in each length bin at each sampling occasion in each gear - unaged
	4darray PltM(1,nlakes,1,meffort,1,ngears,1,nbins);			//number of recaps captured in each length bin at each sampling occasion in each gear - unaged
	5darray PlatU(1,nlakes,1,meffort,irec,ObsAge,1,ngears,1,nbins);			//number of unmarks captured in each length bin at each sampling occasion in each gear - aged
	5darray PlatM(1,nlakes,1,meffort,irec,ObsAge,1,ngears,1,nbins);			//number of recaps captured in each length bin at each sampling occasion in each gear - aged
	
	matrix Eout(1,nlakes,1,nyears+irec);
	3darray Nout(1,nlakes,1,nyears*ndays,irec,ObsAge);
	3darray Lout(1,nlakes,1,Amax,1,nyears);
	3darray L2out(1,nlakes,1,nyears*ndays,irec,ObsAge);
	
	number fpen;
	number decoy;
	LOC_CALCS
		if(dend4!=999)	{cout<<"parameter input error; dend= "<<dend4<<endl;	exit(1);}
		if(debug==1)	cout<<"parameter section done"<<endl;
		cout<<"Initial Calculations Complete"<<endl;
	END_CALCS

PRELIMINARY_CALCS_SECTION
    Gconfill();			if(debug==1)	cout<<"gconfill complete"<<endl;
    LWFitter();		if(debug==1)	cout<<"LW fitter complete"<<endl;

PROCEDURE_SECTION
	fpen=0.;	decoy=0.;
	pM=mfexp(lM);
	pq=mfexp(lq);
	for(int lake=1;lake<=nlakes;lake++)
	{
		//pq(lake)=mfexp(lq);
		palpha(lake)=mfexp(lqalpha);
		pbeta(lake)=mfexp(lqbeta);	//mfexp(lqbeta*10.);
		//pgamma(lake)=1./(1.+mfexp(-tqgamma));
		pgamma(lake,1)=1./(1.+mfexp(-tqgamma(1)));
		pgamma(lake,3)=1./(1.+mfexp(-tqgamma(2)));
		pgamma(lake,2)=1./(1.+mfexp(- -20.));
	}
	pCVl=mfexp(lCVl);
	sigtheta=mfexp(lsigtheta);
	//ptheta=mfexp(ltheta);
	rtheta(1)=Ctheta;
	rtheta(2)=Caththeta;
	rtheta(3)=Dtheta;
	rtheta(4)=Megtheta;
	rtheta(5)=Mtheta;
	rtheta(6)=Ntheta;
	rtheta(7)=Ptheta;
	rtheta(8)=Wtheta;
	int j=1;
	for(int lake=1;lake<=nlakes;lake++)
	{
		int i=1;
		ptheta(lake)=mfexp(mutheta(lake)+rtheta(lake,1)-0.5*square(sigtheta(lake)));
		int st=thetayears(lake)-intheta(lake)+1;
		for(int year=st;year<=thetayears(lake);year++)
		{	ptheta(lake,year)=mfexp(mutheta(lake)+rtheta(lake,i)-0.5*square(sigtheta(lake)));
			if(sd_phase())	logD(j)=log(ptheta(lake,year));
			i++;	j++;
		}
	}

	age_structure();		if(debug==1)	cout<<"age structure complete"<<endl;
	likelihood();				if(debug==1)	cout<<"likelihood complete"<<endl;
//	exit(1);

//_________________________________________
FUNCTION Gconfill
		Gcon.colfill(1,Gconinit);
		for(int lake=1;lake<=nlakes;lake++)
		{
			for(int year=2;year<=tnyears(lake);year++)
			{
				Gcon(lake,year)=Gconinit(lake)*Gconmults(lake,gconind(lake,year-1));
			}
		}

//_________________________________________
FUNCTION LWFitter		//LWFitter function from Carl and is exactly the model used in the previous paper (length-year growth)
	double fs=0.1;
	for(int lake=1;lake<=nlakes;lake++)
	{
		Lage(lake)=0.;
		dvar_matrix Wage(1,Amax,1,nyears(lake));									//modelled weight at age for each lake-year combo
		Lage(lake,1)=LengthHatch;		Wage(1)=alw*pow(LengthHatch,plw);
		for(int year=1;year<=tnyears(lake);year++)
		{
			for(int a=2;a<=int(Amax);a++)
			{
				int nstep=2;
				dvector ty(1,nstep+1);
				double par_dt=DayT/nstep;
				dvariable L,w, ws;
				dvariable fmult, fmultcon, gro, der;

				dvariable Ra, Efac, relalloc, ders, dlast, dlasts, Pg;
				double tyr;
				ty.fill_seqadd(Ti(a-1),par_dt);
				if(ty(1)-int(ty(1))<=DayT && year>1)
				{	L=Lage(lake,a-1,year-1);	w=Wage(a-1,year-1);}
				else
				{	L=Lage(lake,a-1,year);	w=Wage(a-1,year);}
				ws=fs*alw*pow(L,plw);
				for(int ti=1;ti<=nstep;ti++)
				{
					tyr=ty(ti)-int(ty(ti));
					double timeofyear=(tyr*365+1);
					fmult=pow(Qten,(Temp(lake,timeofyear)-Tbar(lake))/10.);
					fmultcon=pow(QtenCon,(Temp(lake,timeofyear)-Tbar(lake))/10.)*mfexp(-Ktemp*(Temp(lake,timeofyear)-Tinact))/(1.+mfexp(-Ktemp*(Temp(lake,timeofyear)-Tinact)));
					gro=fmultcon*Gcon(lake,year)*pow(ws/fs,Gpow);
					der=gro-fmult*Mcon*pow(w,Mpow);
					Ra=(1.-w/WinfFactor);
					if(Ra<1e-5)	{Ra=1e-5;}
					Efac=-rs*(fs*w*pow(Ra,WsPower)-ws);
					if(Efac>100)	{Efac=100;}
					if(Efac<(-100))	{Efac=-100;}
					relalloc=rbase/(1.+mfexp(Efac));
					ders=fs*der*relalloc;
					Pg=0.;
					if(tyr<=TySpawn && (tyr+par_dt)>TySpawn)	{Pg=Pgonad*w/(1.+mfexp(-100.*(w-Wma)/Wma));}
					gro=gro*fs+0.001;
					ders=ders/(1.+mfexp(-5.*ders/gro));
					if(ti==1)
					{
						w+=der*par_dt-Pg;
						ws+=ders*par_dt;
					}
					else
					{
						w+=par_dt/2.*(3.*der-dlast)-Pg;
						ws+=par_dt/2.*(3.*ders-dlasts);
					}
					dlast=der;
					dlasts=ders;
					if(w<1e-20)	{w=1e-20;}
					if(ws<1e-20)	{ws=1e-20;}
					if(w>1e6)	{w=1e6;}
					if(ws>1e6)	{ws=1e6;}
					L=pow(ws/(fs*alw),(1./plw));
				}
				Lage(lake,a,year)=L;
				Wage(a,year)=w;
			}
		}
	}
	if(outLen==1)
	{
		Lout=Lage;
		ofstream lens("Rainbow Lengths.txt");
		lens<<Lout<<endl;
	}

  //_________________________________________
FUNCTION age_structure
	if(debug==1)	cout<<"start age structure"<<endl;
	int check=0;		//set check==1 on both the simulator and estimator to check the abundance time series of marked and unmarked fish from year-class 6
	int outlake=2;			//lake to spit out abundance data - only on converged model
	//detail=0;
	PlatU.initialize();		PlatM.initialize();
	PltU.initialize();		PltM.initialize();
	int JanAge,				year;		//JanAge is the age (in 1/DayT's) of each year-class at Jan 1 of the first sampling year
	dvariable mu;//,			pSurv,		mort;
	dvar_vector dist(1,ngroups),			gtgbin(1,ngroups),			pSurv(1,ngroups),			mort(1,ngroups);
	dvector lengths(1,ngroups);
	dvar_vector fullF(1,ngroups);
	dvar_vector shunt(1,studlen);
	shunt=0.;
	Eout=0.;

	dvar_matrix NaU(1,myclass,1,ngroups);									//numbers of unmarked fish at age for each yearclass in each timestep
	dvar_matrix NaM(1,ObsAge,1,ngroups);									//numbers of marked fish at large in each yearclass and timestep
	dvar_matrix NageMark(1,ObsAge,1,ngroups);		//numbers marked from each ageclass from each growth type group
//	dvar_matrix addmark(1,ObsAge,1,ngroups);
//	dvar_matrix submark(1,ObsAge,1,ngroups);
	dvar_matrix AMort(1,ObsAge,1,ngroups);				//numbers that died in each ageclass from each growth type group
	dvar_matrix fates(1,ngears,1,3);
	dvar_matrix vuln(1,ngears,1,ngroups);							//selectivity of each of the three gears, and a forth which is the product of all three gears and used for 'other' gear
	dvar_matrix F(1,ngears,1,ngroups);

	for(int lake=1;lake<=nlakes;lake++)
	{if(detail==1)	cout<<"lake"<<endl;
		NaU=0.;			NaM=0.;
		NageMark=0.;	AMort=0.;
//		addmark=0.;	submark=0.;
		//	INITIALIZE THE POPULATION
		for(int yclass=1;yclass<=thetayears(lake);yclass++)
		{
			JanAge=max(int(ObsAge-yclass+(1.-TyHatch))*ndays,recage);
			year=max(yclass-ObsAge+startyr(lake)+irec-1-(tstart(lake)-1),1);
			mu=Lage(lake,JanAge,max(year-tstart(lake)+1,1));		//mean length of year-class
			double dmu=value(mu);
			lengths=dmu+0.99*dmu*groups;							//growth type groups extend ~3 standard deviations from the mean length at age
			dist=dnorm(lengths,mu,pCVl(lake)*mu);				//normal distribution of gtg's around mean length of year-class  ## FORMERLY DNORM2 FROM MY
			dist/=sum(dist);
			pSurv=1.;														//proportion surviving to Jan 1 of first year of sampling or year of recruitment
			if(yclass<(ObsAge-irec+1))
				pSurv=mfexp(-pM(lake)*pow(lengths/lr,c)*(JanAge-recage)/ndays);
			NaU(max(yclass,1))+=ptheta(lake,max(yclass,1))*A(lake)*elem_prod(dist,pSurv);		//initial number of unmarked fish is estimated density times lake area distributed over growth-type-groups
			//cout<<yclass<<"\t"<<sum(NaU(max(yclass,1)))<<endl;
		}
		if(outLen==1)
		{
			for(int firstyear=1;firstyear<=irec+startyr(lake)-1;firstyear++)
			{
				for(int yclass=1-firstyear;yclass<=ObsAge-firstyear;yclass++)
				{
					int SpawnAge=(ObsAge-yclass-firstyear+(1.+TySpawn-TyHatch))*ndays;
					//cout<<yclass<<"\t"<<SpawnAge<<endl;
					int iage=ObsAge-yclass+1-firstyear;
					mu=Lage(lake,SpawnAge,1);
					double dmu = value(mu);
					lengths=dmu+0.99*dmu*groups;
					dist=dnorm(lengths,mu,pCVl(lake)*mu);
					dist/=sum(dist);
					dvar_vector tempNa=ptheta(lake,max(yclass,1))*A(lake)*elem_prod(dist,mfexp(-pM(lake)*pow(lengths/lr,c)*SpawnAge/ndays));
					dvar_vector pmat=elem_div(pow(lengths,3.776817),pow(lengths,3.776817)+pow(232.0493,3.776817));
					dvar_vector fec=0.034674*pow(lengths,1.85);
					dvariable pfemale=0.5;			//proportion of spawners that are female
					dvariable pviab=0.9;				//proportion of eggs in ovaries that are viable and fertilized (van Winkle et al. 1998)	
					Eout(lake,startyr(lake)+irec-firstyear)+=sum(elem_prod(elem_prod(pmat,fec)*pfemale*pviab,tempNa));
				}
			}
		}
		if(detail==1)	cout<<"start age"<<endl;
		//	PROCEED THROUGH TIME-STEPS BY REMOVING MORTALITIES AND MOVING MARKS BETWEEN SUB-POPULATIONS
		int yclass,			age;
		int starti=1;
		int trip=0;
		for(int year=startyr(lake);year<=nyears(lake);year++)
		{if(detail==1)	cout<<"year\t"<<year<<endl;
			NageMark=0.;	AMort=0.;
//			addmark=0.;	submark=0.;
			for(int iage=irec+1;iage<=ObsAge;iage++)				//iage is integer age
				AMort(iage)=-NaM(iage-1);				//adds the marked fish back to the unmarked fish at the end of the year
			for(int day=1;day<=ndays;day++)
			{if(detail==1)	cout<<"day\t"<<day<<endl;
				//this part adds observed marks to the marked population and removes observed mortalities from the whole population
				for(int iage=irec;iage<=ObsAge;iage++)
				{
					age=(iage-TyHatch)*ndays+day;
					yclass=ObsAge-iage+year-startyr(lake)+1;
					mu=Lage(lake,age,max(year-tstart(lake)+1,1));
					double dmu = value(mu);
					lengths=dmu+0.99*dmu*groups;							//growth type groups extend ~3 standard deviations from the mean length at age
					mort=1.;
					if(trip==0)	mort=mfexp(-(pM(lake)*pow(lengths/lr,c))/36.);
					NaU(yclass)=elem_prod(posfun(NaU(yclass)-AMort(iage),tiny,decoy),mort);					//unmarked fish with marks and unmarked morts removed
					if(day>1)
						NaM(iage)=elem_prod(posfun(NaM(iage)+NageMark(iage),tiny,decoy),mort);			//marked fish with new marks added and marked morts removed
					else	NaM(iage)=0.;
					//cout<<"NaU\t"<<sum(NaU(age,yclass))<<"\t"<<sum(NaM(iage,day))<<endl;
					if(check==1 && yclass==6)	cout<<sum(NaM(iage))<<"\t"<<sum(NaU(yclass))<<endl;
					if(trip==0 && sd_phase())
					{	shunt((lake-1)*8*ndays+(year-1)*ndays+day)+=sum(NaU(yclass)+NaM(iage));}				//	ESTIMATE DENSITY OF FISH
					if(trip==0 && outLen==1)
					{	Nout(lake,(year-1)*ndays+day,iage)=sum(NaU(yclass)+NaM(iage));
						L2out(lake,(year-1)*ndays+day,iage)=sum(elem_prod(square(lengths),NaU(yclass)+NaM(iage)));
						if(day <=TySpawn*ndays && (day+1)>TySpawn*ndays)
						{
							dvar_vector pmat=elem_div(pow(lengths,3.776817),pow(lengths,3.776817)+pow(232.0493,3.776817));
							dvar_vector fec=0.034674*pow(lengths,1.85);
							dvariable pfemale=0.5;			//proportion of spawners that are female
							dvariable pviab=0.9;				//proportion of eggs in ovaries that are viable and fertilized (van Winkle et al. 1998)	
							Eout(lake,year+irec)+=sum(elem_prod(elem_prod(pmat,fec)*pfemale*pviab,NaU(yclass)+NaM(iage)));
						}
					}
				}
				if(detail==1)	cout<<"start catch "<<lake<<" "<<year<<" "<<day<<endl;
				//This part assigns probabilities that observed marks and mortalities are from different age-classes and growth-type groups
				NageMark=0.;		AMort=0.;
//				addmark=0.;	submark=0.;
				int endi=neffort(lake);
				for(int i=starti;i<=endi;i++)
				{
					if((EYear(lake,i)==year && int(EDate(lake,i)*ndays)==day) || trip==1)
					{
						trip=0;
						if(sum(Catch(lake,i))>=1)
						{	//if the modelled day coincides with a day where any fish were caught (either in gear or the 'other' gear) or if this is part of a multi-day depletion
	//						count(i)=0.;
							fates=0.;
//							dvar_vector temptake(1,6);
							PlatU(lake,i)=0.;	PlatM(lake,i)=0.;	PltU(lake,i)=0.;	PltM(lake,i)=0.;
							dvar_vector Ft(1,ngears);
							for(int gear=1;gear<=ngears;gear++)
							{
								Ft(gear)=pq(lake,gear)*Et(lake,gear,i);
								if(sum(Catch(lake,i,gear)(1,3)+Catch(lake,i,gear,6))>0.)	
								{
									fates(gear,1)=sum(Catch(lake,i,gear)(1,3))/(sum(Catch(lake,i,gear)(1,3))+Catch(lake,i,gear,6));	//unmarked fish removed from NaU (dying or being marked)
									fates(gear,2)=(Catch(lake,i,gear,1))/(sum(Catch(lake,i,gear)(1,3))+Catch(lake,i,gear,6));			//unmarked fish added to NaM (being marked and not dying)
								}
								if(sum(Catch(lake,i,gear)(4,5))>0.)
									fates(gear,3)=Catch(lake,i,gear,4)/sum(Catch(lake,i,gear)(4,5));												//recaped fish dying
							}
//							temptake=colsum(Catch(lake,i));
							//cout<<"i\t"<<i<<fates<<endl;
							for(int iage=irec;iage<=ObsAge;iage++)		//This loop creates matrices of v(l)*N(a)*P(l|a)
							{	if(detail==1)	cout<<"lake\t"<<lake<<"\t"<<"age1\t"<<iage<<endl;
								fullF=0.;	F=0.;
								yclass=ObsAge-iage+year-startyr(lake)+1;
								age=(iage-TyHatch)*ndays+day;
								mu=Lage(lake,age,max(year-tstart(lake)+1,1));
								double dmu = value(mu);
								lengths=dmu+0.99*dmu*groups;
								gtgbin=(lengths+dw-min_l)/delt_l;					//distributes gtg's into the appropriate size-bin
								for(int gear=1;gear<=ngears;gear++)
								{
									vuln(gear)=1./(1.-pgamma(lake,gear))*pow((1.-pgamma(lake,gear))/pgamma(lake,gear),pgamma(lake,gear))*elem_div(mfexp(palpha(lake,gear)*pgamma(lake,gear)*(pbeta(lake,gear)-lengths)),(1.+mfexp(palpha(lake,gear)*(pbeta(lake,gear)-lengths))));
									F(gear)=vuln(gear)*Ft(gear);
									fullF+=F(gear);
								}
								fullF=posfun(fullF,1e-20,decoy);
								for(int gear=1;gear<=ngears;gear++)
								{
									F(gear)=elem_prod(elem_div(F(gear),fullF),1.-mfexp(-fullF));
									for(int gtg=1;gtg<=ngroups;gtg++)
									{
										int gbin=value(gtgbin(gtg));
										if(gbin>=1.&gbin<nbins)
										{
											PlatU(lake,i,iage,gear,gbin)+=NaU(yclass,gtg)*F(gear,gtg);
											PlatM(lake,i,iage,gear,gbin)+=NaM(iage,gtg)*F(gear,gtg);
											PltU(lake,i,gear,gbin)+=NaU(yclass,gtg)*F(gear,gtg);
											PltM(lake,i,gear,gbin)+=NaM(iage,gtg)*F(gear,gtg);
										}
									}
									if(detail==1)	cout<<"calc F"<<endl;
									if(detail==1)	cout<<"move fish"<<endl;
									AMort(iage)+=elem_prod(NaU(yclass),F(gear))*fates(gear,1);			//any fish caught from unmarked population is lost to that population
									NageMark(iage)+=elem_prod(NaU(yclass),F(gear))*fates(gear,2)		//fish caught from the unmarked population in fyke nets and marked are gained by the marked population
										-elem_prod(NaM(iage),F(gear))*fates(gear,3);		//fish caught in the marked population that die are lost to that population
//									addmark(iage)+=elem_prod(NaU(yclass),F(gear))*fates(gear,2);
//									submark(iage)+=elem_prod(NaM(iage),F(gear))*fates(gear,3);
								}
								//cout<<"fates\t"<<fates<<"\t"<<sum(column(Catch(lake,i),5))/sum(Catch(lake,i))<<"\t"<<sum(column(Catch(lake,i),6))/sum(Catch(lake,i))<<endl;
								//count(i)+=elem_prod(NaU(yclass)+NaM(iage),fullF);
							}
//							if(sum(AMort)>0.)	{AMort/=sum(AMort);		AMort*=sum(temptake(1,3));}
//							if(sum(addmark)>0.)	{addmark/=sum(addmark);	addmark*=temptake(1);}
//							if(sum(submark)>0.)	{submark/=sum(submark);	submark*=temptake(4);}
//							NageMark=addmark-submark;
							trip=Eflag(lake,i);
							day-=Eflag(lake,i);
							endi=starti;
							starti=i+1;
							//exit(1);
						}
					}
				}
			}
		}
	}
	if(sd_phase())	
		preddens=log(shunt+1e-10);
	//cout<<preddens<<endl;
	if(check==1)	
		exit(1);
	if(outLen==1)
	{
		ofstream nums("Rainbow Abundance.txt");
		nums<<Nout<<endl;
		ofstream L2s("Rainbow Effective Density.txt");
		L2s<<L2out<<endl;
		ofstream fec("Rainbow Egg Deposition.txt");
		fec<<Eout<<endl;
		cout<<"Lengths, Numbers and Effective Densities output"<<endl;
		exit(1);
	}

//_________________________________________
FUNCTION likelihood
	dvariable Ulf,			Uaf,			Mlf,			Maf;
	Ulf=0.;					Uaf=0.;		Mlf=0.;		Maf=0.;
	//detail=0;
	double aprop=0.;
	dvariable sdpen=0.;
	for(int lake=1;lake<=nlakes;lake++)
	{
		for(int i=1;i<=neffort(lake);i++)
		{
			for(int gear=1;gear<=ngears;gear++)
			{	if(detail==1)	cout<<"i "<<i<<"\t"<<Et(lake,gear,i)<<"\t"<<sum(Catch(lake,i,gear))<<endl;
				if(Et(lake,gear,i)>1e-10 && sum(Catch(lake,i,gear))>=1)
				{
					aprop=1.-sum(CUnmarkFreq(lake,gear,i,mage)+CMarkFreq(lake,gear,i,mage))/sum(CUnmarkFreq(lake,gear,i)+CMarkFreq(lake,gear,i));
					Ulf+=dlpois(CUnmarkFreq(lake,gear,i,mage)/Et(lake,gear,i),posfun(PltU(lake,i,gear)*(1.-aprop)/Et(lake,gear,i),tiny,fpen));//+tiny);
					//cout<<"unmarked\t"<<i<<"\t"<<gear<<"\t"<<"trip "<<Eflag(lake,i)<<"\n"<<CUnmarkFreq(lake,gear,i,mage)<<"\n"<<PltU(lake,i,gear)<<"\n"<<Ulf<<endl;
					Mlf+=dlpois(CMarkFreq(lake,gear,i,mage)/Et(lake,gear,i),posfun(PltM(lake,i,gear)*(1.-aprop)/Et(lake,gear,i),tiny,fpen));//+tiny);
					//cout<<"marked\t"<<i<<"\t"<<gear<<"\t"<<"trip "<<Eflag(lake,i)<<"\n"<<CMarkFreq(lake,gear,i,mage)<<"\n"<<PltM(lake,i,gear)<<"\n"<<Mlf<<endl;
					if(aprop>1e-10)
					{
						for(int iage=irec;iage<=ObsAge;iage++)
						{
							for(int l=min_l;l<=nbins;l++)
							{
								if(CUnmarkFreq(lake,gear,i,iage,l)+CMarkFreq(lake,gear,i,iage,l)>0.5)
								{
									Uaf+=dlpois(CUnmarkFreq(lake,gear,i,iage,l)/Et(lake,gear,i),PlatU(lake,i,iage,gear,l)*aprop/Et(lake,gear,i)+tiny);
									Maf+=dlpois(CMarkFreq(lake,gear,i,iage,l)/Et(lake,gear,i),PlatM(lake,i,iage,gear,l)*aprop/Et(lake,gear,i)+tiny);
								}
							}
						}
					}
				}
			}
		}
		sdpen+=sum(lsigtheta(lake)+square(rtheta(lake))/(2.*square(sigtheta(lake))));
	}if(detail==1)	exit(1);
	//sdpen=sum(lsigtheta+square(rtheta)/(2.*square(mfexp(lsigtheta))));
	//exit(1);
	
	//--------------- PRIORS ON GROWTH PARAMETERS ---------------//
	
	dvar_vector gro_priors(1,ngro_par);
	gro_priors(1)=dnorm(Mcon,gro_mean(1),gro_sd(1));
	gro_priors(2)=full_dnorm(Qten,gro_mean(2),gro_sd(2));
	gro_priors(3)=full_dnorm(QtenCon,gro_mean(3),gro_sd(3));
	for(i in 1:nlakes)
	{
		int j = i+3;
		gro_priors(j)=full_dnorm(Gconinit(i),gro_mean(j),gro_sd(j));
	}
	for(i in 1:gconlen)
	{
		int j = i+3+nlakes;
		gro_priors(j)=full_dnorm(Gconmults(i),gro_mean(j),gro_sd(j));
	}

	//sdpen+=sum(log(thetasd(1))+square(Ctheta)/(2.*square(thetasd(1))));
	dvariable penmult=1.*fpen;		//penalty should probably be turned off (set to zero)
	f-=Uaf+Maf+Ulf+Mlf-penmult-sdpen;
	f-=sum(full_dnorm(pM,0.58,0.058));
	if(est==1)	
		cout<<"Uf "<<Ulf<<" "<<Uaf<<"\t"<<"Mf "<<Mlf<<" "<<Maf<<"\t"<<"fpen "<<penmult<<" "<<sdpen<<" "<<"f "<<f<<endl;
	//cout<<pq<<endl;

//_________________________________________
FUNCTION dvar_vector P_l_a(const dvariable& mu, const dvariable & std)
	{
		dvar_vector A(bins.indexmin(),bins.indexmax());
		for(int i=bins.indexmin();i<=bins.indexmax();i++)
			A[i]=cumd_norm((bins[i]+dw-mu)/std)-cumd_norm((bins[i]-dw-mu)/std);
		return(A/sum(A));
	}

GLOBALS_SECTION
//	#include "MyLikelihoods.cpp"
	#include "statsLib.h"
	#include <admodel.h>
	#include <df1b2fun.h>
	#include <adrndeff.h>
	
REPORT_SECTION

