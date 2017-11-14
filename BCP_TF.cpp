#include	<iostream>
#include	<fstream>
#include	<string>
#include	<vector>
#include	<cmath>
#include	<iomanip>
#include	<stdlib.h>
#include	<stdio.h>
#include	<unistd.h>
#include	"./TNT/tnt.h"
#include	"./MyFun/myFun.h"
#include	"cppoisson_TF.cpp"

using namespace std;
using namespace TNT;


int sum_dif,shift_size,num_data_frag,num_input_frag,num_block,num_seg,num_peak;
double bg_line,cutline;


const int L1=3000000;
const int L2=1500000;
const int N0=1;
const int N=26;

vector<vector<int> > Plus_data(N);
vector<vector<int> > Minus_data(N);
vector<vector<int> > Plus_input(N);
vector<vector<int> > Minus_input(N);

int *data_frag;
int *input_frag;
double (*data)[5];
double *jump;
double *fold;
double (*ss)[6];
double (*pc)[6];
double (*peak)[9];
double st[100001];

void printhelp();
int cmp ( const void *a , const void *b );
void ct(vector<int>data1,int thred1, vector<int> & out1);
void diff(vector<int> loc_plus1, vector<int> loc_minus1, vector<int>& loc_dif1);
void insert(int a, int b, vector<vector<int> > &I);
void seg(double thre,double (*input)[5]);
void peak_candidate(double t,double (*input1)[5],double(*input2)[6],double(*input3));
double choose_max(int);
int frag_count(int);
double est_lam(int,int);
double prob_pois(int,double);

int needhelp =0;
int next_option;
int frag_size = -1;
//int frag_size2 = -1;
int fe=10;
double p_value = 1e-8;
int multiresults = 0;

char * ChIP_data =0,
     * Input_data =0,
	 * Peak_Results =0;



int main(int argc, char* argv[])
{	
	while ( (next_option = getopt(argc, argv, "h:1:2:e:p:3:m:")) != -1) {
        int this_option_optind = optind ? optind : 1;
        switch (next_option) {
		
			case 'h':	sscanf (optarg,"%d",&needhelp);	break;
			case '1':	ChIP_data = optarg;  break;
			case '2':	Input_data = optarg;  break;

			case 'e':	sscanf (optarg,"%d",&fe);		break;
			case 'p':	sscanf (optarg,"%lf",&p_value); break;

			case '3':	Peak_Results = optarg;	break;
			case 'm':	sscanf (optarg,"%d",&multiresults);	break;

			default:        break;
		}
	}

  if(needhelp ==1) printhelp();
  else
  {	
	FILE *out1,*out2,*out3,*out4,*out5,*out6,*out7,*out8,*out9,*out10,*out11,*out12;
	
	int i,j,chr;
	string filename,test_name, chrName, info,c1;
	char chrSign,name4print[30];
	int m,chrNum,a,b,pre,num,t1,t2;
	double bg,f_plus,f_minus,l1,l2;
	double num_temp,temp,sum,len,average,alpha, beta,thre,c,max1,max2,pc_len,lambda,fold_cutline;


/*******Distribute data to each chromosome******/
	filename = ChIP_data;
	ifstream fin(filename.c_str());
	if(!fin.is_open())
	{
		 cout<<"The ChIP-seq data set doesn't exit."<<endl;
        	 exit(1);
	}
	
 
	while(!fin.eof())
	{
		fin>>chrName>>a>>b>>info>>c1>>chrSign;
		if(fin.fail()) break;
		if(frag_size<0) frag_size = b-a;

		test_name = chrName.substr(0,3);
		if (test_name!="chr") {cout<<"The ChIP-seq data is not in a right format."<<endl; exit(1);}
		chrName=chrName.substr(chrName.find('r')+1,chrName.size()-1);
		if(chrName=="x"||chrName=="X") chrNum=23;
		else if (chrName=="y" || chrName=="Y") chrNum=24;
		else if (chrName=="M" || chrName=="m") chrNum=25;
		else chrNum=atoi(chrName.c_str());

		if(chrNum<N)
		{
			if(chrSign=='+') Plus_data[chrNum].push_back(a);
			else if(chrSign=='-') Minus_data[chrNum].push_back(a);
		}
	}
	//cout<<"read size:"<<frag_size<<endl;

	fin.close();
	fin.clear();

	filename=Input_data;
	fin.open(filename.c_str());
	if(!fin.is_open()) {
		cout<<"The control data set  doesn't exit."<<endl; 
		exit(1);
	}

	
	while(!fin.eof())
	{
		fin>>chrName>>a>>b>>info>>c1>>chrSign;
		if(fin.fail()) break;

		test_name = chrName.substr(0,3);
		if (test_name!="chr") {cout<<"The control data is not in a right format."<<endl; exit(1);}

		chrName=chrName.substr(chrName.find('r')+1,chrName.size()-1);
		if(chrName=="x"||chrName=="X") chrNum=23;
		else if (chrName=="y" || chrName=="Y") chrNum=24;
		else if (chrName=="M" || chrName=="m") chrNum=25;
		else chrNum=atoi(chrName.c_str());

		if(chrNum<N)
		{
			if(chrSign=='+') Plus_input[chrNum].push_back(a);
			else if(chrSign=='-') Minus_input[chrNum].push_back(a);
		}
	}
	fin.close();
	fin.clear();

	cout<<"finish distribute the data"<<endl;




/*******Data: estimate shift size, reorder and remove noise ******/	
	vector <int> loc_dif;
	for (chr=N0;chr<N;chr++)
	{
		int * temp1;
		int l=Plus_data[chr].size();
		if(l>0)
		{
			
			temp1 = new int[l];
			for(i=0; i<l; i++)
				temp1[i] = Plus_data[chr][i];
			qsort(temp1,l,sizeof(int),cmp);
			Plus_data[chr].clear();
			pre=0;
			for (i=0; i<l; i++)
			{
				if(temp1[i]!=pre) Plus_data[chr].push_back(temp1[i]);
				pre = temp1[i];
			}
			delete []temp1;
		}

		l=Minus_data[chr].size();
		
		if(l>0)
		{
			temp1 = new int[l];
			for(i=0; i<l; i++)
				temp1[i] = Minus_data[chr][i];
			qsort(temp1,l,sizeof(int),cmp);
			Minus_data[chr].clear();
			pre=0;
			for (i=0; i<l; i++)
			{
				if(temp1[i]!=pre) Minus_data[chr].push_back(temp1[i]);
				pre = temp1[i];
			}	
			delete []temp1;
		}
		if(Plus_data[chr].size()>0&&Minus_data[chr].size()>0)
		{
			bg =(double)Plus_data[chr].size()/(Plus_data[chr][Plus_data[chr].size()-1]-Plus_data[chr][0]);
			f_plus = ceil(fe*1000*bg);

			bg = (double)Minus_data[chr].size()/(Minus_data[chr][Minus_data[chr].size()-1]-Minus_data[chr][0]);
			f_minus = ceil(fe*1000*bg);


			vector <int>loc_plus;
			vector <int>loc_minus;

			ct(Plus_data[chr],f_plus,loc_plus);	



			ct(Minus_data[chr],f_minus,loc_minus);

			diff(loc_plus,loc_minus,loc_dif);
		}
	}
	sum_dif = 0;
	for (i=0;i<loc_dif.size();i++) sum_dif +=loc_dif[i];
	shift_size = sum_dif/(2*loc_dif.size());

	cout<<"shift size:"<<shift_size<<endl;


	out1 = fopen(Peak_Results,"w");
	if (multiresults ==1)
	{        
		out2 = fopen("results_1e-5","w");
        out3 = fopen("results_5e-6","w");
        out4 = fopen("results_e-6","w");
        out5 = fopen("results_5e-7","w");
        out6 = fopen("results_e-7","w");
        out7 = fopen("results_5e-8","w");
        out8 = fopen("results_e-8","w");
        out9 = fopen("results_5e-9","w");
        out10 = fopen("results_e-9","w");
        out11 = fopen("results_5e-10","w");
        out12 = fopen("results_e-10","w");
	}

	for (chr=N0;chr<N;chr++)
	{
		if (Plus_data[chr].size()+Minus_data[chr].size()>0)
		{
			/***********shift data**********/
			data_frag = new int[Plus_data[chr].size()+Minus_data[chr].size()+1] ;	
			data_frag[0]=0;
			if (Plus_data[chr].size()>0)
                for (i=0; i<Plus_data[chr].size(); i++) data_frag[i+1]=Plus_data[chr][i]+shift_size;
			num=Plus_data[chr].size()+1;

			if (Minus_data[chr].size()>0)
			{
				for (i=0; i<Minus_data[chr].size(); i++)
					if(Minus_data[chr][i]>shift_size) {data_frag[num]=Minus_data[chr][i]-shift_size;num++;}
					cout<<"number of data chr"<<chr<<":"<<num-1<<endl; 
			}
			
			num_data_frag = num-1;
			qsort(data_frag,num_data_frag+1,sizeof(int),cmp);

//////////////////////***************************************************///////////////////////////////
			input_frag = new int[Plus_input[chr].size()+Minus_input[chr].size()+1] ;
			input_frag[0]=0;	
			if (Plus_input[chr].size()>0)
			{
				int * temp1;
				temp1 = new int[Plus_input[chr].size()];
				for(i=0;i<Plus_input[chr].size();i++) temp1[i] = Plus_input[chr][i]+shift_size;
				qsort(temp1,Plus_input[chr].size(),sizeof(int),cmp);

				pre=-1;num=1;
				for (i=0; i<Plus_input[chr].size(); i++)
				{
					if(temp1[i]!=pre) {input_frag[num]=temp1[i];num++;}
					pre = temp1[i];
				}			
				delete []temp1;
			}

			if (Minus_input[chr].size()>0)
			{
				int * temp1;
				temp1 = new int[Minus_input[chr].size()];
				for(i=0;i<Minus_input[chr].size();i++) temp1[i]= Minus_input[chr][i]-shift_size;
				qsort(temp1,Minus_input[chr].size(),sizeof(int),cmp);
				pre=-1;
				for (i=0; i<Minus_input[chr].size(); i++)
					if(temp1[i]!=pre && temp1[i]>=0) {input_frag[num]=temp1[i];num++;pre = temp1[i];}
				cout<<"number of input chr"<<chr<<":"<<num-1<<endl; 
				delete []temp1;
			}
			num_input_frag = num-1;// in fact the num_data should be num-1, just for convenient for qsort;
			qsort(input_frag,num_input_frag+1,sizeof(int),cmp);
			cout<<"finish pre-processing data of chr"<<chr<<endl;

			/*************calculate block data,posterior mean,ss,pc and peak****************/
			vector<vector<int> > block;
			for(i=1;i<=num_data_frag;i++)	insert(data_frag[i],data_frag[i]+frag_size,block);
			num_block = block.size();
			cout<<"number of blocks:"<<num_block<<endl;

			Matrix <int> obs(num_block+1,3);
			data = new double [num_block+1][5];
			jump = new double [num_block+1];

			for(i=0; i<num_block;i++)
			{	
				obs[i+1][0]=block[i][0];
				obs[i+1][1]=block[i][1];
				obs[i+1][2]=block[i][2];
			}


			temp = floor(log10((double)num_block)+0.5);
			double p = 1.0/pow(10.0,temp*1.0);

			sum=len=0.0;
			for (i=2;i<=num_block;i++){
				len += (obs[i][1]-obs[i][0]+1.0);
				sum += obs[i][2]*(obs[i][1]-obs[i][0]+1.0);
			}
			average = sum/len;

			if(num_block<7) {cout<<"Too less data for the model."<<endl;exit(1);}
			
			cppoisson  tmp(obs,p,1.0,1.0,7,4);
			tmp.BcmixSmooth();

			i=0;
			data[i][0]=data[i][1]=data[i][2]=data[i][3]=jump[i]=0.0; 
			for(i=1;i<=num_block;i++){
				data[i][0] = obs[i][0];
				data[i][1] = obs[i][1];
				data[i][2] = data[i][1] - data[i][0]+1.0;
				data[i][3] = obs[i][2];
				data[i][4] = tmp.estPara[i];
				jump[i] = data[i][4]-data[i-1][4];
			}

			cout<<"Finish calculate the posterior mean "<<endl;

			st[0]=0;                  //calculating the factorial first
			for (i=1;i<=100000;i++)
				st[i] = st[i-1] + log(i);
			i=0;
			while(prob_pois(i,average)<0.999) i++;
			bg_line = (i>=1)?i:1;


			ss = new double[num_block][6];
			seg(bg_line,data);


			pc =new double [num_seg+1][6];
			peak_candidate(1.0,data,ss,jump);


			peak = new double[num_seg+1][9];
			fold = new double [num_seg+1];fold[0] = 0.0;

			for (i=1;i<=num_seg;i++){
				peak[i][0]=pc[i][0];
				peak[i][1]=pc[i][1];
				peak[i][2]=pc[i][2];
				peak[i][3]=pc[i][3];
				peak[i][4]=peak[i][5]=0;
				peak[i][6]=peak[i][7]=1.0;
				temp = 0.0;
				for (j=pc[i][4];j<=pc[i][5];j++){
					if (temp<=data[j][3]) {temp = data[j][3]; peak[i][8] = data[j][1];}
				}
			}
			
			c=(double)num_data_frag/(double)num_input_frag;
			if(c>10){cout<<"The ChIP-seq data is much more larger than the input data. It might cause unreliable results."<<endl;}

			bg=(num_input_frag-1)/(double)(input_frag[num_input_frag-1]+frag_size-1-input_frag[1]);

			for (i=1;i<=num_seg;i++)
			{
				pc_len = peak[i][2];
				max1 = bg*pc_len;

				l1 = frag_count(i);
				peak[i][4] = (double)l1;

				max2 = choose_max(i)*pc_len;
				l2 = (max1>=max2)?max1:max2;
				peak[i][5] = c*l2;


				peak[i][6] =fabs(1.0-prob_pois(l1,peak[i][5]));

				peak[i][7] = peak[i][4]/peak[i][5];
				fold[i] = peak[i][7]; 

			}

			qsort((void*)fold,num_seg+1,sizeof(double),cmp);		
			if (num_seg<=3000) fold_cutline=2;
			else{
				fold_cutline = fold[num_seg-3000]>=2?fold[num_seg-3000]:2;
				fold_cutline = fold_cutline<=6?fold_cutline:6;
			}
	
			if (chr>=1&&chr<=22) sprintf(name4print,"chr%d",chr);
			if (chr==23) sprintf(name4print,"chrX");
			if (chr==24) sprintf(name4print,"chrY");
			if (chr==25) sprintf(name4print,"chrM");

			for (i=1;i<=num_seg;i++)
				if((peak[i][6]<=p_value)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
					fprintf(out1, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);


			if (multiresults==1)
			{
				for (i=1;i<=num_seg;i++)
				{
					if((peak[i][6]<=1.0e-05)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out2, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=5.0e-06)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out3, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=1.0e-06)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out4, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=5.0e-07)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out5, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=1.0e-07)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out6, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=5.0e-08)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out7, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=1.0e-08)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out8, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=5.0e-09)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out9, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=1.0e-09)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out10, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=5.0e-10)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out11, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);

					if((peak[i][6]<=1.0e-10)&&(peak[i][4]>peak[i][5])&&(peak[i][7]>=fold_cutline))
						fprintf(out12, "%s\t%d\t%d\t%d\t%lf\t%e\t%d\n",name4print,(int)(peak[i][0]-1),(int)peak[i][1],(int)peak[i][2],peak[i][3],peak[i][6],(int)peak[i][8]);
				}

			}
			cout<<"Finish choosing the significant peaks."<<endl;

			delete []data;
			delete []jump;
			delete []ss;
			delete []pc;
			delete []peak;
			delete []fold;
			delete []data_frag;
			delete []input_frag;
		}
		else cout<<"chr"<<chr<<"is empty."<<endl;
	}
  	fclose(out1);

	if(multiresults==1)
	{
		fclose(out2);
        	fclose(out3);
        	fclose(out4);
        	fclose(out5);
        	fclose(out6);
        	fclose(out7);
        	fclose(out8);
        	fclose(out9);
        	fclose(out10);
        	fclose(out11);
      		fclose(out12);	
	}
  }
	return 0;
}

void printhelp()
{
    cout<<"\t-1\tThe ChIP-seq data set you want to input."<<endl;
    cout<<"\t-2\tThe control/input data set you want to input."<<endl;
    cout<<"\t-e\tThe fold enrichment parameter helps to estimate the shift size. Range: 5---12, default is 8."<<endl;
    cout<<"\t-p\tThe p_value you want to use for remove false positive based on control data. Range:1e-5---1e-15, default is 1e-8."<<endl;
    cout<<"\t-3\tThe results data set with 7 columns you want to output."<<endl;
	cout<<"\t-m\tThe flag indicates whether you want to output multiple results corresponding to different p_value. 1 means output the multiple results while 0 means don't output. Default: 0."<<endl;
	cout<<"\t-h\tThe flag indicates whether you want to output the help manual. 1 means outputting the help instructions while 0 means NOT outputting. Default: 0. If you set as 1, the program would not works except displaying the manual."<<endl;
}

int cmp ( const void *a , const void *b )
{
	return *(int *)a - *(int *)b;
}

double prob_pois(int l1, double lambda){
	int i1;
	double prob,logprob,sum_temp;
	prob = exp(-1*lambda);
	for (i1=1;i1<=l1;i1++){
		logprob=0.0;
		sum_temp = st[i1];
		logprob = (-lambda)+i1*log(lambda)-sum_temp;
		prob += exp(logprob);
	}
	return prob;
}

void ct(vector<int>data1,int thred1, vector<int> & out1)
{
	int loop,rec,rec1,rec2,loc,len1,max,loc_max;

	rec1=0;rec2=1;
	while (rec2<data1.size()-1)
	{
		rec = rec1;
		while((data1[rec]+frag_size)>data1[rec2]&&rec2<data1.size()) 
		{rec= rec2;rec2++;}
		
		if ((rec2-rec1)>=thred1)
		{
			len1 = data1[rec]+frag_size-data1[rec1];
	
			vector<int> count(len1+1);
			for (loop=rec1;loop<rec2;loop++)
			{
				for(loc=1;loc<=frag_size;loc++) count[data1[loop]-data1[rec1]+loc]++;
			}
			max = count[0];
			for (loc=0;loc<=len1;loc++) 
			{
				if (max<count[loc])
				{	
					max=count[loc];
					loc_max = loc+data1[rec1];
				}
			}
			if (max >=thred1) out1.push_back(loc_max);
		}
		rec1 = rec2;rec2++;
	}

}

void diff(vector<int> loc_plus1, vector<int> loc_minus1, vector<int> & loc_dif1)
{
    int i=0,j=0;
    int temp_diff;
    bool b=false;

    for(; i<loc_plus1.size()&!b&j<loc_minus1.size();i++)
    {
        while(loc_plus1[i]>=loc_minus1[j]&&j<loc_minus1.size()-1) j++;

	temp_diff = loc_minus1[j]-loc_plus1[i];

       	if(temp_diff<200 && temp_diff>10)
        {
	        loc_dif1.push_back(temp_diff);
		j++;
	}
	if(j>=loc_minus1.size()||temp_diff<=0)   
		b=true;
    }
}

void insert(int a, int b, vector<vector<int> > & I)
{
	vector<int> v(4);
	v[0]=1;
	v[1]=a;
	v[2]=0;
	v[3]=a;
	vector<vector<int> >::iterator it;
	if(I.size()<1)
	{
		I.push_back(v);
	}
	if(I.size()==1) 
	{
		v[0]=a+1;
		v[1]=b;
		v[2]=1;
		v[3]=b-a;
		I.push_back(v);
	}
	else
	{
		it=I.end();
		it--;
		int tail=(*it)[1];
		if(tail<a)
		{
			v[0]=tail+1;
			v[1]=a;
			v[2]=0;
			v[3]=a-tail;
			I.push_back(v);
			v[0]=a+1;
			v[1]=b;
			v[2]=1;
			v[3]=b-a;
			I.push_back(v);
		}
		else if(tail==a)
		{
			if((*it)[2]==1)
			{
				(*it)[1]=b;
				(*it)[3]=b-(*it)[0]+1;
			}
			else 
			{
				v[0]=a+1;
				v[1]=b;
				v[2]=1;
				v[3]=b-a;
				I.push_back(v);
			}
		}
		else//tail>a
		{
			/*add the tail part to the I as one more element*/
			if(tail!=b)
			{
				/*add the tail part to the I as one more element*/
				v[0]=tail+1;
				v[1]=b;
				v[2]=1;
				v[3]=b-tail;
				I.push_back(v);
				it=I.end();
				it=it-2;
			}
			else
			{
				it=I.end();
				it=it-1;
			}
			bool stop=false;//stop is true, when it reachs the begining of I or the change is done; false, overwise.
			while(tail>a&&!stop)
			//stop when the current tail is larger than a ((a,b) will not affect this element) or the begining of I is reached.
			{
				
				int head=(*it)[0];
				//cout<<"head:"<<head<<endl;
				if(head<a+1)
				{
					/*if part of the element is covered by (a,b), then the current element is the first element invovled*/
					/*record the new element between the first and second invovled elements as v*/
					v[0]=a+1;
					v[1]=(*it)[1];
					v[2]=(*it)[2]+1;
					v[3]=(*it)[1]-a;
					/*change the values of the first element invovled*/
					(*it)[1]=a;
					(*it)[3]=a-head+1;
					/*call the installed method of vector:insert*/
					it=I.insert(it+1,v);
					/*move it back to the begining of the first element invovled*/
					it--;
					stop=true;
				}
				else if(head==a+1 && (*it)[2]==(*(it-1))[2]-1)
				{
					(*(it-1))[1]=tail;
					(*(it-1))[3]=tail-(*(it-1))[0]+1;
					it=I.erase(it);
					stop=true;
				}
				else
				{
					/*if the whole element is covered by (a,b), raise the count by 1*/
					(*it)[2]++;
				}
				
			
				if(it!=I.begin())
				{
					/*if it is not the begining of I, move forwards by 1*/
					it--;
					tail=(*it)[1];
				}
				else stop=true;
				
			}			
		}
	}
}
/////smooth out the blocks whose posterior mean is less than 1 and combine the remain blocks//////////
void seg(double thre,double (*input)[5]){

	int i,j,k,j1;
	double temp_sum;

	j=1;i=1;

	while (1){
		temp_sum = 0.0;
		if (j >num_block) break;
		if (input[j][4] <= thre){
			j1=j;
			while ((j1<=num_block)&&(input[j1][4]<=thre)) j1++;
			j = j1 ;
		}
		else {
			j1=j;
			while ((j1<=num_block)&&(input[j1][4] > thre)) j1++;
			ss[i][0] = input[j][0];
			ss[i][1] = input[j1-1][1];
			ss[i][2] = ss[i][1] -ss[i][0] + 1.0;
    			for (k=j; k<j1; k++) 
				temp_sum += input[k][2]*input[k][4];
			ss[i][3] = temp_sum/ss[i][2];
			ss[i][4] = j;
			ss[i][5] = j1-1;
			i++;
			j = j1;

		}
	}
	num_seg = i-1;
}

void peak_candidate(double t,double (*input1)[5],double(*input2)[6],double(*input3)){

	int j,k,j1,j2;
	double max;
	int l1,l2,max_start,max_end;
	
	pc[0][0]=pc[0][1]=pc[0][2]=pc[0][3]=pc[0][4]=pc[0][5]=0.0;
	for (j=1;j<=num_seg;j++){
			
			j1 = (int)input2[j][4];
			j2 = (int)input2[j][5];

			max = input1[j1][4];

			if (j2 == j1){
				pc[j][0] = input1[j1][0];
				pc[j][1] = input1[j1][1];
				pc[j][2] = input1[j1][2];
				pc[j][3] = input1[j1][4];
				pc[j][4] = pc[j][5] = j1;
			}
                        
			else{
				for(k=j1;k<=j2;k++)
					if(max<input1[k][4]) max = input1[k][4];

				l1 = j1; l2 = j2;
				while (input1[l1][4]!=max) l1++;
				max_start = l1;
				while (input1[l2][4]!=max) l2--;
				max_end = l2;
				
				
				while((l1>j1)&&(fabs(input3[l1])<t)) l1--;
				max_start = l1;
				while((l2<j2)&&(fabs(input3[l2+1])<t)) l2++;
				max_end=l2;

				pc[j][0] = input1[max_start][0];
				pc[j][1] = input1[max_end][1];
				pc[j][2] = pc[j][1]-pc[j][0]+1.0;
				pc[j][3] =  max;
				pc[j][4] = max_start;
				pc[j][5] = max_end;

			}
	}
}

double choose_max(int rec){
	int length;
	int start[6];
	int end[6];
	double lam[6];
	double max=0.0;
	length = (int)peak[rec][2];

	for(int i=0;i<6;i++)
	{
		start[i]=(int)peak[rec][0]-i*length;
		end[i]=(int)peak[rec][1]+i*length;
		lam[i]=est_lam(start[i],end[i]);
		if(max<=lam[i]) max=lam[i];
	}
	
	return max;
}
int frag_count(int rec){
        int start_ind,end_ind,comp_ind;
        int ind1,ind2;
        int count;

        start_ind =1;
        end_ind = num_data_frag;
        while (start_ind<end_ind-1){
			comp_ind=start_ind+(end_ind-start_ind)/2;
            if(data_frag[comp_ind]==peak[rec][1]){
				start_ind = comp_ind;
			break;
		}
		else{
			if(data_frag[comp_ind]<peak[rec][1])
				start_ind=comp_ind;
			else
				end_ind=comp_ind;
		}

	}
	ind1 = start_ind;


	start_ind =1;
	end_ind = num_data_frag;
	while(start_ind<end_ind-1){
		comp_ind=start_ind+(end_ind-start_ind)/2;
		if((data_frag[comp_ind]+frag_size) ==peak[rec][0]){
			end_ind = comp_ind;
			break;
		}
		else{
			if((data_frag[comp_ind]+frag_size)>peak[rec][0])
				end_ind=comp_ind;
			else start_ind = comp_ind;
		}
	}
	ind2=end_ind;

	count = ind1-ind2+1;
	return count;
}



double est_lam(int st,int ed){
	int start_ind,end_ind,comp_ind;
	int ind1,ind2;   //the last index of the fragment
	int count,len;
	double lam;
	
	len = ed-st+1;
		

	start_ind =1;
	end_ind = num_input_frag;
	while (start_ind<end_ind-1){
		comp_ind=(start_ind+end_ind)/2;
		if(input_frag[comp_ind]==ed){
			start_ind = comp_ind;
			break;
		}
		else{
			if(input_frag[comp_ind]<ed)
				start_ind=comp_ind;
			else
				end_ind=comp_ind;
		}
	}
	ind1 = start_ind;



	start_ind =1;
	end_ind = num_input_frag;
	while(start_ind<end_ind-1){
		comp_ind=(start_ind+end_ind)/2;
		if((input_frag[comp_ind]+frag_size) ==st){
			end_ind = comp_ind;
			break;
		}
		else{
			if((input_frag[comp_ind]+frag_size)>st)
				end_ind=comp_ind;
			else 
				start_ind = comp_ind;
		}
	}       
	ind2=end_ind;

	count = ind1-ind2+1;

	lam = ((double)count)/((double)len);	
	return lam;
}
