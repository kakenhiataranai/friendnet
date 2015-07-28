//統計量算出

#ifndef _CSTATH_
#define _CSTATH_

//statistics: EMアルゴリズムを追加する
//サロゲートデータ生成とwaylandstatisticsを実装


//rank計算用クラス///////////////////////////////////////////
template<class S> class Ranktmp{
  public:
   S d;
   S *p;
   uint r;
};

//ランク計算用比較演算
template<class S> 
bool sortbyd(Ranktmp<S> i,Ranktmp<S> j) {return (i.d<j.d);}
template<class S> 
bool sortbyp(Ranktmp<S> i,Ranktmp<S> j) {return (i.p<j.p);}
//rank計算用関数　///////////////////////////////////////////

class Statistics{
   //dsfmt_t dsfmt;
   vector<int> onlinehistgramtmp;
   double onlinehistgramtmpdatamin, onlinehistgramtmpdatamax;
   double onlinehistgramtmpbinstep;

 public:
   Statistics();
   Statistics(unsigned int initrandom);

   //実行開始時間
   int startingtime;

   //getting statistics
   double average(const vector<double> *);
   double weightaverage(const vector<double> *, const int);
   double variance(const vector<double> *);
   double sum(const vector<double> *);
   double standarddeviation(const vector<double> *);
   double stdev(const vector<double> *);
   double median(const vector<double> *);
   vector<unsigned int> histgram(const double, const double, const double, 
                        const vector<double> *);
   vector<double> histgramdensity(const double, const double, const double, 
                        const vector<double> *);
   double range(const double, const vector<double> *);
   double range_k(const double, const vector<double> *);
   void probabilitynormalize(vector<double> *);
   double autocorr(const vector<double> *, int);
   vector<vector<double> > makebootstrapdata(vector<double> *);
   double regression_linear(vector<double> *, vector<double> *);
   double corrcoef(vector<double> *, vector<double> *);


   //online histgram
   void onlinehistgramset(const double, const double, const double);
   void onlinehistgramput(const double);
   vector<int> onlinehistgramget();

   //distance
   template<class T> double entropy(vector<T> *);
   template<class T> double kullbackleibler(vector<T> *, vector<T> *);
   double kantrovichmetric(vector<double> *, vector<double> *);

   //random values
   //double randuniform(const double);
   //unsigned int irand();
   double randnormal(const double);
   double randchisq(const unsigned int);
   double randlognormal(const double, const double);
   double randt(const unsigned int);
   double randsincos(void);
   double randline(const double);
   double randsin(const double);
   double randcircle(void);
   double randexp(double);
   double randmaxwell(const double);
   double randgamma(double,double);
   double randgamma_old(double);
   void   randdir (double *, double *, int, double);
   double randcauchy(double, double);
   double randinvnormal(double, double);
   double randtriangular(const double);
   double randpower(double);
   double randbeta(double, double);
   double randlogistic(void);
   double randf(double, double);
   double randweibull(double);
   template<class T> int randvector(const vector<T> *, double);
   int randpoisson(double);
   int randgeometric(double);
   int randbinomial(int, double);
   vector<double> randmcmc(uint,uint,double,double (*)(double));

   //value output
   double valnormal(const double, const double, const double);
   double vallognormal(const double, const double, const double);

   //Fitting
   
   //others
   template<class T> vector<T> cv_old(int, vector<T> *);
   template <class T> vector<vector<T> > cv(vector<T>, double, Statistics *);

   double norm(vector<double> *, vector<double> *, int);

   template<class T> vector<T> freerun(vector<T> *, vector<T> *, int, 
                       T (*model) (vector<T> *, vector<T> *));
   vector<double> lsq(vector<double> *, vector<double> *, int);
   template<class T> vector<uint> rank(vector<T> *);

//inline
   template <class T>
   T max(const vector<T> *in){
      return *max_element(in->begin(), in->end());
   }
   template <class T>
   T min(const vector<T> *in){
      return *min_element(in->begin(), in->end());
   }

};
//////////////////////////////////////////////////////////////////////
Statistics::Statistics(){
   startingtime=time(0);
   //dsfmt_gv_init_gen_rand(startingtime);
   //dsfmt_init_gen_rand(&dsfmt, startingtime);
}

Statistics::Statistics(unsigned int initrandom){
   startingtime=time(0);
   //dsfmt_gv_init_gen_rand(initrandom);
   //dsfmt_init_gen_rand(&dsfmt, initrandom);
}

///////////////////////////////////////////////////////////////
//平均
double Statistics::average(const vector<double> *in){
   double out=0.;
   if(in->size()<=1){cout << "Error at statistics." <<endl;}
   out=(double)accumulate(in->begin(),in->end(),0.0);
   out=out/in->size();
   return out;
}

//Lnノルム
double Statistics::weightaverage(const vector<double> *in, int l){
   double out=0.;
   if(in->size()<=1){cout << "Error at statistics." <<endl;}

   for(int i=0;i<in->size();i++){out+=pow((*in)[i],(double)l);}
   out=pow(out, 1./(double)l);
   out=out/in->size();
   return out;
};

//分散
double Statistics::variance(const vector<double> *in){
   double out=0.;
   if(in->size()<=1){cout << "Error at statistics." <<endl;}
   
   out=average(in);
   vector<double> sqr=(*in);
   for(vector<double>::iterator i=sqr.begin();i!=sqr.end();i++){
      (*i)=((*i)-out)*((*i)-out);
   }
   out=(double)accumulate(sqr.begin(),sqr.end(),0.0)/(in->size()-1);
   return out;
}

//和
double Statistics::sum(const vector<double> *in){
   double out=0.;
   if(in->size()<=1){cout << "Error at statistics." <<endl;}
   out=(double)accumulate(in->begin(),in->end(),0.0);
   return out;
}

//標準偏差
double Statistics::standarddeviation(const vector<double> *in){
   double out=0.;
   if(in->size()<=1){cout << "Error at statistics." <<endl;}
   out=variance(in);
   out=sqrt(out);
   return out;
}

//標準偏差略称
double Statistics::stdev(const vector<double> *in){
   return standarddeviation(in);
}


//メジアン
double Statistics::median(const vector<double> *in){
   double out;
   vector<double> tmp=*in;
   sort(tmp.begin(),tmp.end());
   if(tmp.size()%2){out=tmp[tmp.size()/2];}
   else{
      out=(tmp[tmp.size()/2-1]+tmp[tmp.size()/2])/2.;
   }
   return out;
}

//ヒストグラム:単純カウント datamin=datamaxなら範囲自動設定
vector<unsigned int> Statistics::histgram(
                        const double datamin, const double datamax, 
                        const double binstep, const vector<double> *in){
   vector<unsigned int> out;
   unsigned int binsize;
   double datamintmp=datamin, datamaxtmp=datamax;

   if(datamintmp==datamaxtmp){
      datamintmp=min(in);
      datamaxtmp=max(in);
   }

   //Scott's method
   if(binstep==0.){binsize=(int)(3.5*stdev(in)*pow(in->size(),(-1./3.)));}
   else{binsize=(unsigned int)((datamaxtmp-datamintmp)/binstep)+1;}

   out.resize(binsize);
   for(unsigned int i=0;i<binsize;i++){out[i]=0;}
   int size=in->size();
   for(int i=0;i<size;i++){
      if((*in)[i]>=datamintmp && (*in)[i]<=datamaxtmp){
         out[int(((*in)[i]-datamintmp)/binstep)]++;
      }
   }
   return out;
}

//ヒストグラム:密度分布
vector<double> Statistics::histgramdensity(const double datamin, 
                           const double datamax, 
                           const double binstep, const vector<double> *in){

   vector<unsigned int> tmp=histgram(datamin, datamax, binstep, in);
   vector<double> out(tmp.size());
   for(int i=0;i<tmp.size();i++){out[i]=(double)tmp[i]/binstep;}
   return out;
}

double Statistics::range(const double cover, const vector<double> *in){
//カバーする範囲のxの半値幅を返す : cover=0.95なら 例えば-2～2→2を返す

   int halfsize=int(cover*(double)in->size()/2.);
   vector<double> tmp=*in;
   sort(tmp.begin(),tmp.end());
   return (tmp[tmp.size()/2+halfsize]-tmp[tmp.size()/2-halfsize])/2.;
}

double Statistics::range_k(const double k, const vector<double> *in){
//kで示された範囲のxの半値幅を返す
   return k*standarddeviation(in);
}

void Statistics::probabilitynormalize(vector<double> *in){
   //double minval=min(in);

   //for(int i=0;i<in->size();i++){(*in)[i]-=minval;}
   double sum=accumulate(in->begin(),in->end(),0.);
   for(vector<double>::iterator i = in->begin(); i != in->end(); i++ ){
   (*i)/=sum;
   }
}

//自己相関 遅れn
double Statistics::autocorr(const vector<double> *in, int n){
   double out=0.;
   unsigned int size=in->size()-n;

   double mu=average(in);
   double sigma2=variance(in);

   for(int i=0;i<size;i++){
   out+=(((*in)[i]-mu)*((*in)[i+n]-mu)/sigma2);
   }
   return out/size;
}

//ブートストラップのために1点抜いたデータを多数作成
vector<vector<double> > Statistics::makebootstrapdata(vector<double> *in){
   vector<vector<double> >out;
   out.resize(in->size());
   vector<double> mabiki;
   mabiki.resize(in->size()-1);
   int cnt;

   for(int i=0;i<in->size();i++){//outloop
      cnt=0;
      for(int j=0;j<in->size();j++){
         if(i!=j){mabiki[cnt]=(*in)[j];cnt++;}
      }
      out[i]=mabiki;
   }
   return out;
}

//回帰直線の傾き: 底上げ値はyの平均値
double Statistics::regression_linear(vector<double> *x, vector<double> *y){
   double avx=average(x);
   double avy=average(y);

   double tmp1=0., tmp2=0.;
   for(int i=0;i<x->size();i++){
      tmp2+=((*x)[i]-avx)*((*x)[i]-avx);
      tmp1+=((*x)[i]-avx)*((*y)[i]-avy);
   }

   return tmp1/tmp2;
}

//相関係数
double Statistics::corrcoef(vector<double> *x, vector<double> *y){
   uint n=x->size();
   int i;
   double sx, sy, sxx, syy, sxy, dx, dy;

   sx = sy = sxx = syy = sxy = 0;
   for (i = 0; i < n; i++) {
      sx += (*x)[i];  sy += (*y)[i];
   }
   sx /= n;  sy /= n;
   for (i = 0; i < n; i++) {
      dx = (*x)[i] - sx;  dy = (*y)[i] - sy;
      sxx += dx * dx;  syy += dy * dy;  
      sxy += dx * dy;
   }
   sxx = sqrt(sxx / (n - 1));
   syy = sqrt(syy / (n - 1));
   sxy /= (n - 1) * sxx * syy;
   return sxy;
}




//////////////////////////////////////////////////////////////////////////
//非メモリ蓄積型ヒストグラム

// set
void Statistics::onlinehistgramset(const double datamin, const double datamax, 
                                   const double binstep){
   onlinehistgramtmpdatamin=datamin;
   onlinehistgramtmpdatamax=datamax;
   onlinehistgramtmpbinstep=binstep;
   unsigned int binsize=(unsigned int)((datamax-datamin)/binstep)+1;
   onlinehistgramtmp.resize(binsize, 0);
};

// put
void Statistics::onlinehistgramput(const double in){
   if(in>=onlinehistgramtmpdatamin && in<=onlinehistgramtmpdatamax){
      onlinehistgramtmp[
         int((in-onlinehistgramtmpdatamin)/onlinehistgramtmpbinstep)
                       ]++;
   }
};

// get and clear
vector<int> Statistics::onlinehistgramget(){
   vector<int> out=onlinehistgramtmp;
   onlinehistgramtmp.clear();
   return out;
};

//////////////////////////////////////////////////////////////////////////

//エントロピー
template<class T> double Statistics::entropy(vector<T> *in){
   double out=0.;
   T zero=0;
   double sumdiv=1.0/(double)accumulate(in->begin(),in->end(),zero);
   for(vector<double>::iterator i = in->begin(); i != in->end(); i++ ){
      out+=(-(*i)*sumdiv*log((*i)*sumdiv));
   }
   return out;
}

template<class T> double Statistics::kullbackleibler(vector<T> *in1,
                                                     vector<T> *in2){
   double out=0.;
   unsigned int size=in1->size();
   for(int cnt=0;cnt<size;cnt++){
      if((*in1)[cnt]==0.){(*in1)[cnt]+=1e-10;}
      if((*in2)[cnt]==0.){(*in2)[cnt]+=1e-10;}
      out+=((*in1)[cnt]*log((double)(*in1)[cnt]/(double)(*in2)[cnt]));
   }
   return out;
}

//同じxをもつ確率変数のカントロビッチメトリックを求める
double Statistics::kantrovichmetric(vector<double> *in1, vector<double> *in2){
   unsigned int size=in1->size();
   if(size != in2->size()){
      cout << "Error at km. Different size of inputs" <<endl;exit(-1);
   }
   double out=0.;

   vector<double> tmp1, tmp2;
   tmp1.resize(size);
   tmp2.resize(size);
   partial_sum(in1->begin(), in1->end(), tmp1.begin());
   partial_sum(in2->begin(), in2->end(), tmp2.begin());

   for(int i=0; i<size; i++){out+=fabs(tmp1[i]-tmp2[i]);}
   return out;
}


/////////////////////////////////////////////////////////////
//各種乱数
//double Statistics::randuniform(const double i){
//   return fabs(i)*dsfmt_genrand_open_open(&dsfmt);
//   //return fabs(i)*dsfmt_gv_genrand_open_open();
//}

//unsigned int Statistics::irand(){
//   return dsfmt_genrand_uint32(&dsfmt);
//}

double Statistics::randnormal(const double sigma){//E=0, V=sigma
   double rand1=randuniform(1), rand2=randuniform(1);
   return (sqrt(-2.*log(rand1))*sin(2.*PI*rand2)*sigma);
}

double Statistics::randchisq(const unsigned int df){
   double out=0;
   double tmp;
   for(int i=0;i<df;i++){
      tmp=randnormal(1.);
      out+=(tmp*tmp);
   }
   return out;
}

double Statistics::randlognormal(const double mu, const double sd){
   return exp(randnormal(sd)+mu);
};

double Statistics::randt(const unsigned int df){
   return randnormal(1.)/sqrt(randchisq(df)/(double)df);
}

double Statistics::randsincos(void){
   double x=0.0,y=2.0;
   while(y>(2.0*cos(x)*sin(x))){x=randuniform(PI/2.);y=randuniform(2.0);}
   return x;
}

//randline y=2/in/in x
double Statistics::randline(const double in){
   double x=0,y=10;
   double indiv=1./in;
   while(y>(2.0*x*indiv*indiv)){
      x=randuniform((double)in);y=randuniform(2.0*indiv);
   }
   return x;
}

//randsin
double Statistics::randsin(const double in){
   double x=0.0,y=2.0;
   while(y>sin(x)){x=randuniform(in);y=randuniform(1.0);}
   return x;
}

//randcircle
double Statistics::randcircle(void){//-1～1までの円形確率を返す
   double x=0.0,y=2.0;
   while(y>sqrt(1.-x*x)){x=randuniform(2.)-1.;y=randuniform(1.0);}
   return x;
}

//指数分布
double Statistics::randexp(double lambda){
   return -log(randuniform(1.))/lambda;
}

//randmaxwell
double Statistics::randmaxwell(const double T){
   double mdiv=1./4.65e-23;
   double m=4.65e-23;
   double k=1.38e-23;
   double v=0.0,y=1.0;
   double f=0;
   //v_average=sqrt(2kT/m)

   while(y>f){
      v=randuniform(5.*sqrt(2.*k*T*mdiv));y=randuniform(1.0);
      f=4.0*PI*v*v
        *pow(m/(2.0*PI*k*T),1.5)
        *exp((-m*v*v)/(2.0*k*T));
   }
   return v;
}

//randgamma: ガンマ分布
double Statistics::randgamma(double theta, double kappa ){
    int int_kappa;
    double frac_kappa;

    int_kappa  = (int)kappa;
    frac_kappa = kappa - (double)int_kappa;
    
    double u,uu;
    double b,p,x_frac,x_int;
    int i;
    
    //Int
    x_int=0;
    for(i=0;i<int_kappa;i++){
        x_int+=randexp(1.);
    }
    
    //Frac
    if( fabs(frac_kappa) < 0.01 ) x_frac=0;
    
    else{
        b=(exp(1.0)+frac_kappa)/exp(1.0);
        while(1){
            u=randuniform(1.);
            p=b*u;
            
            uu=randuniform(1.);
            
            if(p<=1.0){
                x_frac=pow(p,1.0/frac_kappa);
                if(uu<=exp(-x_frac)) break;
            }
            
            else{
                x_frac=-log((b-p)/frac_kappa);
                if(uu<=pow(x_frac,frac_kappa-1.0)) break;
            }
            
        }
    }
    
    return (x_int+x_frac)*theta;
}

//randgamma_old: ガンマ分布 （旧） mean=1 fixed
double Statistics::randgamma_old(double a){
   double x, y, z;
   double u, v, w, b, c, e;
   int accept = 0;
   if (a < 1){
      /* a < 1. Johnk's generator. Devroye (1986) p.418 */
      e = randexp(1.);
      do {
         x = pow(randuniform(1.), 1 / a);
         y = pow(randuniform(1.), 1 / (1 - a));
      } while (x + y > 1);
      return (e * x / (x + y));
   } else {
      /* a >= 1. Best's rejection algorithm. Devroye (1986) p.410 */
      b = a - 1;
      c = 3 * a - 0.75;
      do {
         /* generate */
         u = randuniform(1.);
         v = randuniform(1.);
         w = u * (1 - u);
         y = sqrt(c / w) * (u - 0.5);
         x = b + y;
         if (x >= 0)
         {
            z = 64 * w * w * w * v * v;
            if (z <= 1 - (2 * y * y) / x)
            {
               accept = 1;
            } else {
               if (log(z) < 2 * (b * log(x / b) - y))
                  accept = 1;
            }
         }
      } while (accept != 1);
      return x;
   }
}

/*
ディリクレ分布からサンプルする場合は, 上の dirrand を使って, double *theta, *alpha をアロケートしてから,
dirrand(theta, alpha, k, 0);
とすれば theta に Dir(alpha) からのサンプルが得られます。最後の引数はディリクレ分布を中心 α/Σ(α) と精度 Σ(α) に分解してサンプリングする時に, 精度(中心への集中の度合い)を調整できるようにするためのパラメータ. 
*/

//ディリクレ分布
void Statistics::randdir (double *theta, double *alpha, int k, double prec){
   double z = 0;
   int i;
   /* theta must have been allocated */
   for (i = 0; i < k; i++)
      if (prec != 0)
         theta[i] = randgamma_old(alpha[i] * prec);
      else
         theta[i] = randgamma_old(alpha[i]);
   for (i = 0; i < k; i++)
      z += theta[i];
   for (i = 0; i < k; i++)
      theta[i] /= z;
}

//コーシー分布
double Statistics::randcauchy(double mu, double gamma){
    return mu + gamma*tan(PI*( randuniform(1.)-0.5 ));
};

//逆正規分布
double Statistics::randinvnormal(double mu, double lambda ){
    double y,w,z;
    y=randchisq(1);
    w= mu+0.5*y*mu*mu/lambda -(0.5*mu/lambda)*sqrt(4.0*mu*lambda*y+mu*mu*y*y); 
    z=randuniform(1.);
    if( z< mu/(mu+w) )  return w;
    else                return mu*mu/w;
}

//三角分布
double Statistics::randtriangular(const double mul){
   return (randuniform(1.)-randuniform(1.))*mul;
}

//累乗分布
double Statistics::randpower(double n){
   return pow(randuniform(1.), n);
}

//β分布
double Statistics::randbeta(double a, double b){
   double x, y;
   do {
      x = pow(randuniform(1.), 1 / a);  y = pow(randuniform(1.), 1 / b);
   } while (x + y > 1);
   return x / (x + y);
}

//logistic分布
double Statistics::randlogistic(void){
   double r;
   r = randuniform(1.);
   return log(r / (1 - r));
}

//F分布
double Statistics::randf(double n1, double n2){
   return (randchisq(n1) * n2) / (randchisq(n2) * n1);
}

//ワイブル分布
double Statistics::randweibull(double alpha){
   return pow(-log(1 - randuniform(1.)), 1 / alpha);
}

//randvector: vectorで与えられた関数に従って乱数発生
template<class T> 
int Statistics::randvector(const vector<T> *in,double maxin){
   //maxinが0でないならそれを優先する
   double max;
   if(maxin==0.0){max=*max_element(in->begin(), in->end())*1.1;}
   else{max=maxin;}
   
   int x=0;
   double y=max;
   
   while(y>(double)((*in)[x])){
      x=(int)randuniform((double)in->size());
      y=randuniform(max);
   }
   return x;
}

//Poisson (ポアソン) 分布
int Statistics::randpoisson(double lambda){
   int k;
   lambda = exp(lambda) * randuniform(1.);
   k = 0;
   while (lambda > 1) {
   lambda *= randuniform(1.);  k++;
   }
   return k;
}

//幾何分布
int Statistics::randgeometric(double p){
   return ceil(log(1 - randuniform(1.)) / log(1 - p));
}

//二項分布
int Statistics::randbinomial(int n, double p){
   int i, r;
   r = 0;
   for (i = 0; i < n; i++)
      if (randuniform(1.) < p) r++;
   return r;
}

//MCMCを使った乱数系列の出力
vector<double> Statistics::randmcmc(uint n, uint burn, double sigma, double (*f)(double)){
Statistics s;
   vector<double> data(n-burn);
   double x=randuniform(1.)-0.5, xnew;

   for(uint i=0;i<n;i++){
      xnew=s.randnormal(sigma)+x;
      if(f(xnew)/(f(xnew)+f(x))>randuniform(1.)){x=xnew;}
      if(i>burn){data[i-burn]=x;}
   }
   return data;
}



////////////////////////////////////////////////////////////////////
//値の出力

double Statistics::valnormal(
 const double mu, double sigma, const double x){
   return 1./sqrt(2.*PI)/sigma*exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
};

double Statistics::vallognormal(
 const double mu, double sigma, const double x){
   return 1./sqrt(2.*PI)/sigma/x
          *exp((-(log(x)-mu)*(log(x)-mu)/(2.*sigma*sigma)));
};


////////////////////////////////////////////////////////////////////

//クロスバリデーションのために配列の一部n個をランダムに切り出す
template<class T> vector<T> Statistics::cv_old(int n, vector<T> *in){
   if(in->size()<n){n=in->size;}
   int skip=(int)randuniform((double)(in->size()-n));
   vector<T> out;
   out.resize(n);
   for(int i=skip;i<n;i++){out[i]=(*in)[skip+i];}
   return out;
}

//クロスバリデーション
template <class T>
vector<vector<T> > Statistics::cv(vector<T> in, double rate, Statistics *s){
   vector<T> out1((uint)(in.size()*rate)), out2(in.size()-out1.size());
   //cout << out1.size() << ": " << out2.size() <<endl;
   //random_shuffle(in.begin(),in.end(), s->randuniform(1.));
   random_shuffle(in.begin(),in.end());
   for(uint i=0;i<out1.size();i++){out1[i]=in[i];}
   for(uint i=0;i<out2.size();i++){out2[i]=in[i+out1.size()];}
   vector<vector<T> > out(2);
   out[0]=out1;
   out[1]=out2;
   return out;
}

//L_nノルムを計算する
double Statistics::norm(vector<double> *in1, vector<double> *in2, int n){
   if(in1->size()!=in2->size()){cout << "Error at norm." <<endl;exit(-1);}
   if(n<1){cout << "Error at norm." <<endl;exit(-1);}

   int size=in1->size();
   double l=0;
   for(int i=0;i<size;i++){
      l+=pow((*in1)[i]-(*in2)[i],n);
   }
   return pow(l, 1./n);
};

//フリーラン
template<class T>
vector<T> Statistics::freerun(vector<T> *coff, vector<T> *init, 
                       int len, 
                       T (*model) (vector<T> *, vector<T> *)){
   vector<T> out(len);
   vector<T> inittmp=(*init);
   
   for(int i=0;i<len;i++){
      out[i]=(*model)(coff,&inittmp);
      //if(out[i]>9e99){out[i]=9e99;}
      rotate(inittmp.begin(), inittmp.begin()+1, inittmp.end());
      inittmp[inittmp.size()-1]=out[i];
   }
   return out;
}

//以下freerunのサンプル
/*
#include "C.h"
#include "C-statistics.h"
double logistic(vector<double> *coef, vector<double> *init){
   double x=(*init)[0];
   return (*coef)[0] *x*(1.-x);
}
int main(){
   C c;
   Statistics s;
   vector<double> init(1);
   init[0]=0.1;
   vector<double> coef(1);
   coef[0]=3.8;
   vector<double> kekka=s.freerun(&coef, &init, 30, &logistic);
   c.vectoroutv(cout, &kekka);
   return 0;
}
*/

//n次最小二乗法(吐き出し法利用)：xはソートされていること
//データ列x,y, 多項式の次数Nを与える: 定数項から高次の順に係数を返す
vector<double> Statistics::lsq(vector<double> *x, vector<double> *y, int N){
   unsigned int n=N+1;
   unsigned int S=x->size();
   unsigned int T=y->size();
   if(S != T){cout << "Error at lsq" << endl; exit(-1);}
   vector<double> out(n);
   vector<vector<double> >A(n);
   for(int i=0;i<n;i++){A[i].resize(n+1,0.0);}

   //初期行列作成
   for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {
         for(int k=0;k<S;k++) {
            A[i][j]+=pow((*x)[k],i+j);
         }
      }
   }

   for(int i=0;i<n;i++) {
      for(int k=0;k<S;k++) {
         A[i][n]+=pow((*x)[k],i)*(*y)[k];
      }
   }

   //吐き出し法
   int pivot;
   double p,q,m,b[1][n+1];

   for(int i=0;i<n;i++){
      m=0;
      pivot=i;
      for(int l=i;l<n;l++) {
         if(fabs(A[l][i])>m) {   //i列の中で一番値が大きい行を選ぶ
            m=fabs(A[l][i]);
            pivot=l;
         }
      }
      if(pivot!=i) {           //pivotがiと違えば、行の入れ替え
         for(int j=0;j<n+1;j++) {
            b[0][j]=A[i][j];   
            A[i][j]=A[pivot][j];
            A[pivot][j]=b[0][j];
         }
      }
   }
   for(int k=0;k<n;k++) {
      p=A[k][k];         //対角要素を保存
      A[k][k]=1;         //対角要素は１になることがわかっているから
      for(int j=k+1;j<n+1;j++) {
         A[k][j]/=p;
      }
      for(int i=k+1;i<n;i++) {
         q=A[i][k];
         for(int j=k+1;j<n+1;j++) {
            A[i][j]-=q*A[k][j];
         }
      A[i][k]=0;         //０となることがわかっているところ
      }
   }

   //解の計算
   for(int i=n-1;i>=0;i--) {
      out[i]=A[i][n];
      for(int j=n-1;j>i;j--) {
         out[i]-=A[i][j]*out[j];
      }
   }
   return out;
}

//順位を出す: vector<uint> kekka=rank(&in);
template<class T> vector<uint> rank(vector<T> *in){
   vector<Ranktmp<T> > tmp(in->size());
   for(uint i=0;i<in->size();i++){
      tmp[i].d=(*in)[i];
      tmp[i].p=&((*in)[i]);
   }
   sort(tmp.begin(), tmp.end(),sortbyd<T>);
   for(uint i=0;i<in->size();i++){
      tmp[i].r=i;
   }
   sort(tmp.begin(), tmp.end(),sortbyp<T>);
   vector<uint> out(in->size());
   for(uint i=0;i<in->size();i++){
      out[i]=tmp[i].r;
   }
   return out;
}

//////////////////////////////////////////////////////////

#endif

