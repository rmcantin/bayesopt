

////////////////////////////////////////////////////////////////////////////////
class kernel
{
public:
  virtual void setISO_Scale( double theta );
  virtual void setARD_Scales( vectord params );

  virtual double operator()( const vectord &x1, const vectord &x2,
			     size_t grad_index) = 0.0;
  virtual ~kernels() {};
};


////////////////////////////////////////////////////////////////////////////////
class MaternIso1: public kernels
{
public:
  void setISO_Scale( double theta ) {mTheta = theta};

  double operator()( const vectord &x1, const vectord &x2,
		     size_t grad_index )
  {
    double r = norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (param_index == 0) 
      return er;
    else 
      return r*er;
  } 
protected:
  double mTheta;
};


////////////////////////////////////////////////////////////////////////////////
class MaternIso3: public kernels
{
public:
  void setISO_Scale( double theta ) {mTheta = theta};

  double operator()( const vectord &x1, const vectord &x2,
		     size_t grad_index )
  {
    double r = sqrt(3) * norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (param_index == 0) 
      return (1+r)*er;
    else 
      return r*r*er; 
  } 
protected:
  double mTheta;
};


////////////////////////////////////////////////////////////////////////////////
class MaternIso5: public kernels
{
public:
  void setISO_Scale( double theta ) {mTheta = theta};

  double operator()( const vectord &x1, const vectord &x2,
		     size_t grad_index )
  {
    double r = sqrt(5) * norm_2(x1-x2)/mTheta;
    double er = exp(-r);

    if (param_index == 0) 
      return (1+r*(1+r/3))*er;
    else 
      return r*(1+r)/3*r*er; 
  } 
protected:
  double mTheta;
};



////////////////////////////////////////////////////////////////////////////////
class SEiso: public kernels
{
public:
  void setISO_Scale( double theta ) {mTheta = theta};

  double operator()( const vectord &x1, const vectord &x2,
		     size_t grad_index )
  {
    double rl = norm_2(x1-x2)/mTheta;
    double k = rl*rl;

    if (param_index == 0)
      return exp(-k/2);
    else 
      return exp(-k/2)*k;
  } 
protected:
  double mTheta;
};


////////////////////////////////////////////////////////////////////////////////
class SEard: public kernels
{
public:
  void setARD_Scales( vectord theta ) {mTheta = theta};

  double operator()( const vectord &x1, const vectord &x2,
		     size_t grad_index )
  {
    vector<double> xd = x1-x2;
    vector<double> ri = ublas_elementwise_div(xd, theta);

    double rl = norm_2(ri);
    double k = rl*rl;

    if (param_index == 0)
      return exp(-k/2);
    else 
      return exp(-k/2)*sqrt(ri(param_index));
  } 
protected:
  vectord mTheta;
};

