#include<iostream>
#include<pthread.h>
#include<cstdlib>

using namespace std;

int count=0;
pthread_mutex_t mutex;

void * threadFcn(void * arg)
{
	int stop=0;

	while (!stop)
	{
		pthread_mutex_lock(&mutex);	
		if (count<1000)
			count++;
		else stop=1;	
		cout << "******** count = " << count << " **********" << endl;
		pthread_mutex_unlock(&mutex);
		if (!stop)
			system("./overdominance");
	}

	return NULL;
}

int main()
{
	const int nthread=7;
	pthread_t threads[nthread];

	pthread_mutex_init(&mutex,NULL);

	for (int i=0;i<nthread;i++)
		int rc=pthread_create(&threads[i],NULL,threadFcn,NULL);

	for (int i=0;i<nthread;i++)
		int rc=pthread_join(threads[i],NULL);

	pthread_mutex_destroy(&mutex);
}

