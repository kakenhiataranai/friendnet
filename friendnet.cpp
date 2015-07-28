//コンストラクタは使わない
//コメントは日本語
//privateは使わない

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cfloat>
#include <climits>
#include <cassert>
#include <cstring>
#include <tr1/random>

#define uint unsigned int

using namespace std;

//グローバルレベル乱数生成
std::tr1::mt19937 mt(time(0));
std::tr1::uniform_real <double> uni;
std::tr1::variate_generator
<std::tr1::mt19937, std::tr1::uniform_real <double> > gen(mt, uni);

//乱数基本関数
double randuniform(double in){
   return in*gen();
}

#include "C-stat.h"

template<class T> void vectorout(ostream& out, const T *in){
   for(uint i=0; i<in->size(); i++){out << (*in)[i] << ",";}
   out <<endl;
}

template<class T> void swap(T *in1, T *in2){
   T tmp=*in1;
   *in1=*in2;
   *in2=tmp;
}

class Person {
  public:
   uint id;
   vector<double> ability;
   vector<double> abilityoriginal;
   set<Person *> friends;
   double effort;
};

class Data{
  public:
   Statistics *s;
   uint iteration;
   uint n;
   uint nettype;
   double netrate;
   double init_ability;
   uint netchange;
   uint N;
   uint outputtype;
   uint timeupdatetype;
   set<pair<Person *, Person *> > links;
   vector<Person> persons;
   map<pair<Person *, Person *>, double> distances;
};


//ランダムネット
set<pair<Person *, Person *> > makenet_random(Data *d){
   uint N=d->persons.size();
   uint linksnum=(uint)(d->netrate*(d->N-1)*d->N/2);

   set<pair<uint, uint> > already;
   set<pair<Person *, Person *> >out;
   uint p1, p2;

   while(1){
      p1=(uint)randuniform(d->N);
      p2=(uint)randuniform(d->N);
      if(p2==p1){continue;}
      if(p2>p1){swap(&p1, &p2);}
      if(already.count(pair<uint, uint>(p1,p2))==1){continue;}
      ((d->persons))[p1].friends.insert(&(((d->persons))[p2]));
      ((d->persons))[p2].friends.insert(&(((d->persons))[p1]));
      already.insert(pair<uint, uint>(p1,p2));
      if(already.size()>=linksnum){break;}
   }

   for(auto i=already.begin();i!=already.end();i++){
      out.insert(pair<Person *, Person *>(&(d->persons[(*i).first]), &(d->persons[(*i).first])));
   }
   return out;
}

//擬似BAモデル
set<pair<Person *, Person *> > makenet_ba(Data *d){
   set<pair<uint, uint> > already;
   set<pair<Person *, Person *> >out;
   uint N=d->persons.size();
   uint linksnum=(uint)(d->netrate*(d->N-1)*d->N/2);
   uint eachlinks;
   vector<uint> counter(1,0);
   uint eachcounter, p;
   
   for(uint i=1;i<d->N;i++){
      eachlinks=(linksnum-already.size())/(d->N-i);//これから作るべきリンク数
      //cout << eachlinks <<endl;
      counter.push_back(0);
      //すべて張る場合
      if(eachlinks>=i){
         for(uint j=0;j<i;j++){
            (d->persons)[j].friends.insert(&((d->persons)[i]));
            (d->persons)[i].friends.insert(&((d->persons)[j]));
            already.insert(pair<uint, uint>(j,i));
            counter[i]++;
            counter[j]++;
         }
      }
      //ランダムに張る場合
      else{
         eachcounter=0;
         while(1){
            p=(uint)d->s->randvector(&counter, 0);
            if(i==p){continue;}
            if(already.count(pair<uint, uint>(p,i))==1){continue;}
            ((d->persons))[i].friends.insert(&(((d->persons))[p]));
            ((d->persons))[p].friends.insert(&(((d->persons))[i]));
            already.insert(pair<uint, uint>(p,i));
            counter[i]++;
            counter[p]++;
            eachcounter++;
            if(eachcounter>=eachlinks){break;}
         }
      }
   }
   for(auto i=already.begin();i!=already.end();i++){
      out.insert(pair<Person *, Person *>(&(d->persons[(*i).first]), &(d->persons[(*i).first])));
   }
   return out;
}

//calc_distances: Person同士の能力差による距離：能力合計は同じである前提
double calc_distances(Data *d, Person *in1, Person *in2){
   uint n=d->n;
   double out=0.;
   for(uint i=0;i<n;i++){
      out+=pow(in1->ability[i]-in2->ability[i],2.);
   }
   return out;
}

//timeupdateのためのサブルーチン
//友人たちの中で最大の能力
double calc_maxfreiendability(Person *p, int n, Data *d){
   double maxval=-0.001;
   for(auto i=p->friends.begin();i!=p->friends.end();++i){
      if((*i)->ability[n]>maxval){maxval=(*i)->ability[n];};
   }
   return maxval;
};

// time_update
void time_update_normal(Data *d){
   vector<double> maxvals(d->n);
   double maxval;
   uint maxloc;
   for(auto i=d->persons.begin();i!=d->persons.end();i++){//各個人に対して

      //友人たちのジャンルごとの能力値最大分布配列を作って自分の能力を引く
      for(uint j=0;j<d->n;j++){//各ジャンルについて
         maxvals[j]=(*i).ability[j]-calc_maxfreiendability(&(*i),j ,d);
      }
      
      //プラスになっているところがなければ、自分の能力に応じて努力を傾ける
      if(d->s->max(&maxvals)<0.){
         maxvals=(*i).ability;
      }
      //プラスになっているところに、その量に応じて努力を配分
      vector<double> abilityplus(d->n, 0.);
      for(uint j=0;j<d->n;j++){
         if(maxvals[j]>0.){abilityplus[j]=maxvals[j];}
      }
      double abilityplussum=d->s->sum(&abilityplus);
      for(uint j=0;j<d->n;j++){
         (*i).ability[j]+=(*i).effort/abilityplussum*abilityplus[j];
      }
   }
   //友人同士の能力の距離を測っておく
   for(auto i=d->persons.begin();i!=d->persons.end();i++){//各個人に対して
      for(auto j=d->persons.begin();j!=i;j++){//各個人に対して
	 d->distances[pair<Person *, Person *>(&(*j),&(*i))]
	   =calc_distances(d, &(*i), &(*j));
      }
   }
}

//他者と能力が差が少ない部分に努力を注ぐ
void time_update_reverse(Data *d){
   vector<double> maxvals(d->n);
   double maxval;
   uint maxloc;
   for(auto i=d->persons.begin();i!=d->persons.end();i++){//各個人に対して

      //友人たちのジャンルごとの能力値最大分布配列を作って自分の能力を引く
      for(uint j=0;j<d->n;j++){//各ジャンルについて
         maxvals[j]=fabs((*i).ability[j]-calc_maxfreiendability(&(*i),j ,d));
      }

      //プラスになっているところに、その量に応じて努力を配分
      vector<double> abilityplus(d->n, 0.);
      for(uint j=0;j<d->n;j++){
         //if(maxvals[j]>0.){
            abilityplus[j]=1./maxvals[j];
         //}
      }
      double abilityplussum=d->s->sum(&abilityplus);
      for(uint j=0;j<d->n;j++){
         (*i).ability[j]+=(*i).effort/abilityplussum*abilityplus[j];
      }
   }
   //友人同士の能力の距離を測っておく
   for(auto i=d->persons.begin();i!=d->persons.end();i++){//各個人に対して
      for(auto j=d->persons.begin();j!=i;j++){//各個人に対して
	 d->distances[pair<Person *, Person *>(&(*j),&(*i))]
	   =calc_distances(d, &(*i), &(*j));
      }
   }
}



void time_reset(Data *d){
    for(auto i=d->persons.begin();i!=d->persons.end();i++){
       (*i).ability=(*i).abilityoriginal;
    }
}

void change_network(Data *d, uint changes){
   vector<pair<Person *, Person *> >dellink;
   vector<pair<Person *, Person *> >addlink;

   for(auto i=d->distances.end();i!=d->distances.begin();i--){
      //もしリンクが既にあったら削除候補にpush_back
      if(d->links.count((*i).first)==1){
         dellink.push_back((*i).first);
      }
      if(dellink.size()>=changes){break;}
   }

   for(auto i=d->distances.begin();i!=d->distances.end();i++){
      //もしリンクがなかったら候補にpush_back
      if(d->links.count((*i).first)==0){
         addlink.push_back((*i).first);
      }
      if(addlink.size()>=changes){break;}
   }
   //実際の切断と結合
   for(auto i=dellink.begin();i!=dellink.end();i++){
      d->links.erase(*i);
      (*i).first->friends.erase((*i).second);
      (*i).second->friends.erase((*i).first);
   }
   
   for(auto i=addlink.begin();i!=addlink.end();i++){
      d->links.insert(*i);
      (*i).first->friends.insert((*i).second);
      (*i).second->friends.insert((*i).first);
   }
}

vector<uint> check_ba(Data *d, uint maxmax){
   vector<double> moto(d->persons.size(), 0.);
   for(uint i=0;i<moto.size();i++){
      moto[i]=(double)d->persons[i].friends.size();
   }
   return d->s->histgram(0.,(double)maxmax,1.,&moto);
}

vector<double> summary_best(Data *d){
   vector<double> out(d->n,0.);
   double maxval;
   //nジャンルそれぞれにおいて、ジャンル中1番の者がもつ能力値（規格化してない）
   for(uint i=0;i<d->n;i++){
      maxval=0.;
      for(uint j=0;j<d->N;j++){
         if(maxval< d->persons[j].ability[i]){maxval=d->persons[j].ability[i];}
      }
      out[i]=(maxval-d->init_ability)/d->iteration;
   }
   return out;
}

vector<double> summary_kl(Data *d, uint histgrams){
   vector<double> out(d->n,0.);
   vector<double> klbase(histgrams, 1./histgrams);
   vector<double> tmp(d->N, 1./d->N);
   vector<double> kl(histgrams);
   for(uint i=0;i<d->n;i++){
      for(uint j=0;j<d->N;j++){
         tmp[j]=(double)d->persons[j].ability[i]/(double)d->iteration;
      }
      kl=d->s->histgramdensity(0.0000000001,1.,1./histgrams, &tmp);
      d->s->probabilitynormalize(&kl);
      out[i]=d->s->kullbackleibler(&kl,&klbase);
   }
   return out;
}

void out_cout(Data *d, vector<double> *s, uint count){
   cout << count << ",";
   cout << d->N << ",";
   cout << d->n << ",";
   cout << d->nettype << ",";
   cout << d->netrate << ",";
   cout << d->iteration << ",";
   cout << d->netchange << ",";
   cout << d->init_ability << ",";
   cout << d->timeupdatetype << ",";
   cout << d->links.size() << ",,";
   d->c->vectorout(cout, s);
}

void out_graphviz(string fn, Data *d){
   ofstream file((fn+".dot").c_str());
   vector<Person> *p=&(d->persons);
   set<pair<Person *,Person *> > *l=&(d->links);
   
   file << "graph G {"<<endl;
   //file << "ranksep=5;\nratio=auto;\n";
   //file << "size = \"10., 10.\";" <<endl;
   for(auto i=p->begin();i!=p->end();i++){
      file << (*i).id <<endl;
   }
   for(auto i=l->begin();i!=l->end();i++){
      file << (*i).first->id <<"--"<< (*i).second->id <<endl;
   }
   file << "}"<<endl;
}

//////////////////////////////////////////////////////
int main(int ac, char *av[]){
   Data d;
   Statistics s;
   d.s=&s;

   if(ac<7){cout 
      << "Op: N n nettype netrate iteration netchange outputtype initability timeupdatetype"
      << endl
      << "Ex: ./ninchi15 300 5 1 0.2 1000 10 1 0.01 1"
      << endl;
      exit(-1);
   }

   d.N=atoi(av[1]);           //
   d.n=atoi(av[2]);           //2-10【変化要素】
   d.nettype=atoi(av[3]);     //【変化要素】
   d.netrate=atof(av[4]);     //0.01-0.9【変化要素】
   d.iteration=atoi(av[5]);   //100-100000
   d.netchange=atoi(av[6]);  //張り替え量
   d.outputtype=atoi(av[7]); //KLか最大値出力か
   d.init_ability=atof(av[8]);//0.01-0.2【変化要素】
   d.timeupdatetype=atoi(av[9]);//normal0かreverse1か

   
   //Person
   d.persons.resize(d.N);
   for(uint i=0;i<d.N;i++){
      d.persons[i].id=i;
      d.persons[i].effort=1.;//【変化要素】
      d.persons[i].ability.resize(d.n); 
      //init ability
      for(uint j=0;j<d.n;j++){
         d.persons[i].ability[j]=randuniform(d.init_ability);
      }
      d.persons[i].abilityoriginal=d.persons[i].ability;
   }

   //make network
   if(d.nettype==0){
      d.links=makenet_random(&d);
   }
   else if(d.nettype==1){
      d.links=makenet_ba(&d);
   }

   //update
   vector<double> summaryresults;
   for(uint ite=0;ite<d.iteration;ite++){
      
      if(d.timeupdatetype==0){time_update_normal(&d);}else
      if(d.timeupdatetype==1){time_update_reverse(&d);}
      
      if(d.outputtype==0){
         summaryresults=summary_best(&d);
      }
      else if(d.outputtype==1){
         summaryresults=summary_kl(&d, 100);
      }
      if(d.netchange!=0){
         change_network(&d, d.netchange);
      }
      out_cout(&d, &summaryresults, ite);
   }
   return 0;
}
