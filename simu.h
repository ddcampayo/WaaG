#ifndef _SIMU_H_
#define _SIMU_H_

class sim_data
{
 public:

  //  sim_data  {}

  int set_no_of_particles(int nn)  {
    h_=std::sqrt(1.0/FT(nn));
    no_of_particles_=nn;
    return nn;
  }

  void set_totalV (const FT totalV ) {
    totalV_ = totalV;
  }

  void set_innerV (const FT innerV ) {
    innerV_ = innerV;
    meanV_  = innerV / FT( no_of_particles_ );
  }

  
  FT totalV() const {
    return totalV_;
  }

  FT meanV() const {
    return meanV_;
  }

  int current_step() const {return current_step_;}

  int next_step()  {return ++current_step_;}

  FT set_dt(FT tt)   {return  dt_ = tt ;}

  FT dt() const  {return  dt_ ;}

  FT advance_time()  {return  time_+=dt() ;}

  FT time() const  {return  time_ ;}

  bool perturb() const {return perturb_;}

  FT do_perturb( FT pp=0 )  { perturb_ = true ;   return pert_rel_ = pp ; }

  FT pert_rel() const {
		    if(perturb_) return pert_rel_;
		    else return 0.0;
  }

  FT h() const { return h_; }

  int no_of_particles() const { return no_of_particles_;}

private:
  int no_of_particles_;
  FT h_;

  int Nsteps_;
  int every_;
  int current_step_;

  FT dt_;
  FT time_;
  
  bool perturb_;
  FT pert_rel_;

  FT totalV_;
  FT innerV_;
  FT  meanV_;

  
};


extern sim_data simu;

#endif
