#include<iostream>
#include<fstream>
#include<boost/random/mersenne_twister.hpp>
#include<boost/random.hpp>

#include<vector>
#include<cmath>

using namespace std;

boost::mt19937 gen(getpid());

struct indiv_t
{
	char gender;
	int pop,gen,mp,fp,id;
	int l1a1,l1a2,l2a1,l2a2;
};

void pick(int &val,int vals[],double probs[])
{
//	double rnd=1.0*rand()/RAND_MAX;
	boost::uniform_01<> dist;

	double rnd=dist(gen);

	double cdf;
	int i=0;
	cdf=probs[i];
	while (rnd>cdf)
	{
		i++;
		cdf+=probs[i];
	}
	val=vals[i];
}

double uniform()
{
        boost::uniform_01<> dist;
	return dist(gen);
}

void print(vector<indiv_t> list)
{
	for (int i=0;i<list.size();i++)
	{
		cout << i << " " << list[i].gender << " ";
		cout << list[i].pop << " " << list[i].gen << " ";
		cout << list[i].id << " " << list[i].mp << " " << list[i].fp;
		cout << " " << list[i].l1a1 << " " << list[i].l1a2 << " ";
		cout << list[i].l2a1 << " " << list[i].l2a2 << endl;
	}
}

void append(vector<indiv_t> &dest,vector<indiv_t> src)
{
	for (int i=0;i<src.size();i++)
		dest.push_back(src[i]);
}

int pick(int n)
{
	return rand() % n;
}

void print(vector<vector<int> > data)
{
	for (int i=0;i<data.size();i++)
	{
		for (int j=0;j<data[i].size();j++)
			cout << data[i][j] << "\t";
		cout << endl;
	}
}

int main()
{
	const int popsize=100,noff=5,ngen=200,nloci=2,reps=1;
	double mort_A1=0,mort_A2=0,recomb_rate=0;
	int male_recomb=0;

	int A1=1,A2=2,B1=10,B2=20;

	int alle_names[nloci][2]={A1,A2,B1,B2};
	double alle_freq[nloci][2]={0.5,0.5,0.2,0.8};

	vector<indiv_t> starting_pop,next_gen;

	starting_pop.resize(popsize);
	for (int i=0;i<popsize;i++)
	{
		if (i<popsize/2)
			starting_pop[i].gender='M';
		else
			starting_pop[i].gender='F';
		starting_pop[i].pop=1;
		starting_pop[i].gen=0;
		starting_pop[i].mp=0;
		starting_pop[i].fp=0;
		starting_pop[i].id=i;
		pick(starting_pop[i].l1a1,alle_names[0],alle_freq[0]);
                pick(starting_pop[i].l2a1,alle_names[0],alle_freq[0]);
                pick(starting_pop[i].l1a2,alle_names[1],alle_freq[1]);
                pick(starting_pop[i].l2a2,alle_names[1],alle_freq[1]);
	}

	print(starting_pop);

	append(next_gen,starting_pop);

	vector<int> gen_tot,het_tot_A,het_tot_B;
	vector<double> alle_freq_tot_A1,alle_freq_tot_A2;

	int g=0;
	while (1)
	{
		g++;
		cout << "g = " << g << endl;
		vector<indiv_t> male_parents,female_parents;
		for (int i=0;i<starting_pop.size();i++)
			if (starting_pop[i].gender=='M')
				male_parents.push_back(starting_pop[i]);
			else
				female_parents.push_back(starting_pop[i]);
		vector<vector<int> > parent_matrix;
		for (int i=0;i<starting_pop.size();i++)
		{
			int mp=pick(male_parents.size());
			int fp=pick(female_parents.size());
			vector<int> tmp_parent;
			tmp_parent.push_back(mp);
			tmp_parent.push_back(fp);
			parent_matrix.push_back(tmp_parent);
		}
		vector<indiv_t> offspring_matrix;
		char gender[]={'M','F'};
		for (int i=0;i<popsize;i++)
		{
			for (int j=0;j<noff;j++)
			{
				int mp=parent_matrix[i][0];
				int fp=parent_matrix[i][1];
				indiv_t offspring;
				offspring.gender=gender[pick(2)];
				offspring.pop=1;
				offspring.gen=g;
				offspring.mp=mp;
				offspring.fp=fp;
				
				int mchrome1[]={starting_pop[mp].l1a1,
						starting_pop[mp].l1a2};
				int mchrome2[]={starting_pop[mp].l2a1,
						starting_pop[mp].l2a2};

				double randad=uniform();
				int mctmp1[2],mctmp2[2];
				if (randad<recomb_rate && male_recomb)
				{
					mctmp1[0]=mchrome1[0];
					mctmp1[1]=mchrome2[1];
					mctmp2[0]=mchrome2[0];
					mctmp2[1]=mchrome1[1];
				}
				else
				{
					mctmp1[0]=mchrome1[0];
					mctmp1[1]=mchrome1[1];
					mctmp2[0]=mchrome2[0];
					mctmp2[1]=mchrome2[1];
				}
				if (pick(2)==0)
				{
					offspring.l1a1=mctmp1[0];
					offspring.l1a2=mctmp1[1];
				}
				else
				{
					offspring.l1a1=mctmp2[0];
					offspring.l1a2=mctmp2[1];
				}

                                int fchrome1[]={starting_pop[fp].l1a1,
                                                starting_pop[fp].l1a2};
                                int fchrome2[]={starting_pop[fp].l2a1,
                                                starting_pop[fp].l2a2};

				double ranmom=uniform();
                                int fctmp1[2],fctmp2[2];
                                if (ranmom<recomb_rate)
                                {
                                        fctmp1[0]=fchrome1[0];
                                        fctmp1[1]=fchrome2[1];
                                        fctmp2[0]=fchrome2[0];
                                        fctmp2[1]=fchrome1[1];
                                }
                                else
                                {
                                        fctmp1[0]=fchrome1[0];
                                        fctmp1[1]=fchrome1[1];
                                        fctmp2[0]=fchrome2[0];
                                        fctmp2[1]=fchrome2[1];
                                }	
				if (pick(2)==0)
				{
					offspring.l2a1=fctmp1[0];
					offspring.l2a2=fctmp1[1];
				}
				else
				{
					offspring.l2a1=fctmp1[0];
					offspring.l2a2=fctmp1[1];
				}
				double rnd=uniform();
				if (offspring.l1a1==A1&&offspring.l2a1==A1
						&&rnd<mort_A1)
					break;
				if (offspring.l1a1==A2&&offspring.l2a1==A2
						&&rnd<mort_A2)
					break;
				offspring.id=offspring_matrix.size();
				offspring_matrix.push_back(offspring);
			}
		}
		vector<indiv_t> maleoff,femaleoff;
		for (int i=0;i<offspring_matrix.size();i++)
			if (offspring_matrix[i].gender=='M')
				maleoff.push_back(offspring_matrix[i]);
			else
				femaleoff.push_back(offspring_matrix[i]);

		starting_pop.clear();
		for (int i=0;i<popsize/2;i++)
		{
			starting_pop.push_back(maleoff[pick(maleoff.size())]);
			starting_pop.push_back(femaleoff[pick(femaleoff.size())]);
		}
		append(next_gen,starting_pop);
		int count_A1=0,count_A2=0,count_B1=0,count_B2=0;
		for (int i=0;i<popsize;i++)
		{
		count_A1+=(starting_pop[i].l1a1==A1)+(starting_pop[i].l2a1==A1);
		count_A2+=(starting_pop[i].l1a1==A2)+(starting_pop[i].l2a1==A2);
		count_B1+=(starting_pop[i].l1a2==B1)+(starting_pop[i].l2a2==B1);
		count_B2+=(starting_pop[i].l1a2==B2)+(starting_pop[i].l2a2==B2);
		}
		double alle_freq_A1=count_A1/(2.0*popsize);
		double alle_freq_A2=count_A2/(2.0*popsize);
		double alle_freq_B1=count_B1/(2.0*popsize);
		double alle_freq_B2=count_B2/(2.0*popsize);
		double het_A=1.0-(pow(alle_freq_A1,2)+pow(alle_freq_A2,2));
		double het_B=1.0-(pow(alle_freq_B1,2)+pow(alle_freq_B2,2));
		het_tot_A.push_back(het_A);
		het_tot_B.push_back(het_B);
		gen_tot.push_back(g);
		alle_freq_tot_A1.push_back(alle_freq_A1);
		alle_freq_tot_A2.push_back(alle_freq_A2);

		cout << alle_freq_A1 << " " << alle_freq_A2 << endl;

/*
		if (alle_freq_A1==0)
		{
			cout << "No more A1 after generation " << g << endl;
			break;
		}
		if (alle_freq_A2==0)
		{
			cout << "No more A2 after generation " << g << endl;
			break;
		}
*/
		if (g>ngen)
			break;
	}
		
}

