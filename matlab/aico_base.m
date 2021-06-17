double sAICO::step(){
  uint t, T=sys->get_T();
  
  rememberOldState();
  
  switch(sweepMode){
    //NOTE: the dependence on (sweep?..:..) could perhaps be replaced by (dampingReference.N?..:..)
    //updateTimeStep(uint t, bool updateFwd, bool updateBwd, uint maxRelocationIterations, bool forceRelocation){
    case smForwardly:
#if 1
      for(t=1; t<=T; t++) updateTimeStep(t, true, false, (sweep?0:1), (cost<0.));
      for(t=T+1; t--;)    updateTimeStep(t, false, true, 1, (true));
#else
      for(t=1; t<=T; t++) updateFwdMessage(t);
      if(cost<0.) for(t=0; t<=T; t++) updateTaskMessage(t, b[t]());
      for(t=T+1; t--;) if(!(fixFinalState && t==T)) updateBwdMessage(t);
      for(t=0; t<=T; t++) updateBelief(t);
      for(t=0; t<=T; t++) updateTaskMessage(t, b[t]()); //compute task message at reference!
      for(t=0; t<=T; t++) updateBelief(t);
#endif
      break;
    case smSymmetric:
      for(t=1; t<=T; t++) updateTimeStep(t, true, false, 1, (cost<0.));
      for(t=T+1; t--;)    updateTimeStep(t, false, true, 1, (true));
      break;
    case smLocalGaussNewton:
      for(t=1; t<=T; t++) updateTimeStep(t, true, false, (sweep?5:1), (cost<0.)); //relocate iteratively on
      for(t=T+1; t--;)    updateTimeStep(t, false, true, (sweep?5:1), (true)); //...fwd & bwd sweep
      break;
    case smLocalGaussNewtonDamped:
      for(t=1; t<=T; t++) updateTimeStepGaussNewton(t, true, false, (sweep?5:1)); //GaussNewton in fwd & bwd sweep
      for(t=T+1; t--;)    updateTimeStep(t, false, true, (sweep?5:1), (true));
      break;
    case smILQG:
      for(t=T+1; t--;) if(!(fixFinalState && t==T)) updateBwdMessage(t);
      for(t=1; t<=T; t++) unrollXhatFwd(t);
      for(t=0; t<=T; t++) updateTaskMessage(t, xhat[t]()); //compute task message at reference!
      break;
    default: HALT("non-existing sweep mode");
  }
  
  b_step=maxDiff(b_old, b);
  dampingReference=b;
  //if(!sys->isKinematic()) soc::getPositionTrajectory(q, b); else q=b;
  
  //for(t=0; t<=T; t++) updateTaskMessage(t, b[t], 1e-8, 1.); //relocate once on fwd & bwd sweep
  
  if(sys->os){
    *sys->os <<"AICO(" <<sys->get_T() <<") " <<std::setw(3) <<sweep <<" time " <<MT::timerRead(false) <<" setq " <<countSetq <<" diff " <<b_step <<" damp " <<damping <<std::flush;
  }

  //if(cost_old>0 && cost>cost_old)
  //cost_old = sys->analyzeTrajectory(b_old, display>0); //this routine calles the simulator again for each time step
  if(!advanceBeliefBeyondXhat || sweepMode==smILQG){
    cost = evaluateTrajectory(b, display>0); //this routine takes the current R, r matrices to compute costs
    //double exact_cost = sys->analyzeTrajectory(b, display>0); //this routine calles the simulator again for each time step
    //sys->costChecks(b);
    //cout <<"DIFF=" <<fabs(cost-exact_cost) <<endl;
  }else{
    cost = analyzeTrajectory(*sys, b, display>0, &cout); //this routine calles the simulator again for each time step
  }
  
  //-- analyze whether to reject the step and increase damping (to guarantee convergence)
  if(damping) perhapsUndoStep();
  
  sweep++;
  displayCurrentSolution();
  return b_step;
}

void AICO::iterate_to_convergence(){
  self->iterations_till_convergence=0;
  for(uint k=0; k<self->max_iterations; k++){
    double d=self->step();
    if(k && d<self->tolerance){
      self->iterations_till_convergence=k+1;
      break; //d*(1.+damping)
    }
  }
}

void sAICO::init_messages(){
  uint T=sys->get_T();
  arr x0;
  sys->get_x0(x0);
  uint n=x0.N;
  //messages
  s.resize(T+1, n);  Sinv.resize(T+1, n, n);
  v.resize(T+1, n);  Vinv.resize(T+1, n, n);
  r.resize(T+1, n);  R.resize(T+1, n, n);     r.setZero();  R   .setZero();
  for(uint t=0;t<=T;t++){
    s[t]=x0;  Sinv[t].setDiag(1e-10);
    v[t]=x0;  Vinv[t].setDiag(1e-10);
  }
  s[0]=x0;  Sinv[0].setDiag(1e10);
  if(useBwdMsg){ v[T] = bwdMsg_v;  Vinv[T] = bwdMsg_Vinv; }

  //beliefs
  b.resize(T+1, n);  Binv.resize(T+1, n, n);  b.setZero();  Binv.setZero();  b[0]=x0;  Binv[0].setDiag(1e10);
  rhat.resize(T+1);     rhat.setZero();
  xhat.resize(T+1, n);  xhat.setZero();  xhat[0]=x0;
  //if(!sys->isKinematic()) soc::getPositionTrajectory(q, b); else q=b;
  //dampingReference = qhat;
  dampingReference.clear();
  
  //resize system matrices
  A.resize(T+1, n, n);  tA.resize(T+1, n, n);  Ainv.resize(T+1, n, n);  invtA.resize(T+1, n, n);
  a.resize(T+1, n);
  if(!sys->isKinematic()){
    B.resize(T+1, n, n/2);  tB.resize(T+1, n/2, n); //fwd dynamics
    Hinv.resize(T+1, n/2, n/2);
  }else{
    B.resize(T+1, n, n);  tB.resize(T+1, n, n); //fwd dynamics
    Hinv.resize(T+1, n, n);
  }
  Q.resize(T+1, n, n);
  
  //initialize system matrices at time 0
  sys->setx(x0);
  sys->getControlCosts(NoArr, Hinv[0](), 0);
  sys->getDynamics(A[0](), tA[0](), Ainv[0](), invtA[0](), a[0](), B[0](), tB[0](), Q[0](), 0);
  
  //delete all task cost terms
  //for(uint i=0;i<phiBar.N;i++){ listDelete(phiBar(i));  listDelete(JBar(i));  }
  phiBar.resize(T+1);
  JBar  .resize(T+1);
  Psi   .resize(T+1, n);  Psi.setZero();
  
  rememberOldState();
}

void sAICO::init(ControlledSystem& _sys){
  sys = &_sys;
  
  MT::getParameter(sweepMode, "aico_sweepMode");
  MT::getParameter(max_iterations, "aico_max_iterations");
  MT::getParameter(maxStepSize, "aico_maxStepSize");
  MT::getParameter(tolerance, "aico_tolerance");
  MT::getParameter(display, "aico_display");
  MT::getParameter(damping, "aico_damping");
  MT::getParameter(advanceBeliefBeyondXhat,"aico_advanceBeliefBeyondXhat");
  
  if(MT::checkParameter<MT::String>("aico_filename")){
    MT::getParameter(filename, "aico_filename");
    cout <<"** output filename = '" <<filename <<"'" <<endl;
    os=new std::ofstream(filename);
  }else{
    os = &cout;
  }
  
  sweep=0;
  scale=0;
  cost=-1;
  useBwdMsg=false;
  fixFinalState=false;
  init_messages();
}

void sAICO::updateTimeStep(uint t, bool updateFwd, bool updateBwd, uint maxRelocationIterations, bool forceRelocation){
  uint T=sys->get_T();
  if(updateFwd) updateFwdMessage(t);
  if(updateBwd) if(!(fixFinalState && t==T)) updateBwdMessage(t);

  updateBelief(t);
  
  for(uint k=0; k<maxRelocationIterations; k++){
    if(k || !forceRelocation)
      if(maxDiff(b[t], xhat[t])<tolerance) break;
    
    updateTaskMessage(t, b[t]());
    
    //optional reupdate fwd or bwd message (if the dynamics might have changed...)
    //if(updateFwd) updateFwdMessage(t);
    //if(updateBwd) updateBwdMessage(t);
    
    if(advanceBeliefBeyondXhat)
      updateBelief(t);
  }
}

void sAICO::updateFwdMessage(uint t){
  CHECK(t>0, "don't update fwd for first time step");
  arr barS, St;
  if(!sys->isKinematic()){
#ifndef TightMode
    inverse_SymPosDef(barS, Sinv[t-1] + R[t-1]);
    St = Q[t-1];
    St += B[t-1]*Hinv[t-1]*tB[t-1];
    St += A[t-1]*barS*tA[t-1];//cout <<endl <<endl <<t <<endl;
    s[t] = a[t-1] + A[t-1]*(barS*(Sinv[t-1]*s[t-1] + r[t-1]));
    inverse_SymPosDef(Sinv[t](), St);
    //cout <<"s\n" <<s[t] <<endl <<Sinv[t] <<endl;
#else
    St = Q[t-1];
    St += B[t-1]*Hinv[t-1]*tB[t-1];
    s[t] = a[t-1] + A[t-1]*qhat[t-1];
    inverse_SymPosDef(Sinv[t](), St);
#endif
  }else{
    inverse_SymPosDef(barS, Sinv[t-1] + R[t-1]);
    s[t] = barS * (Sinv[t-1]*s[t-1] + r[t-1]);
    St = Hinv[t-1] + barS;
    inverse_SymPosDef(Sinv[t](), St);
  }
}

void sAICO::updateBwdMessage(uint t){
  uint T=sys->get_T();
  if(fixFinalState){ CHECK(t!=T, "don't update bwd for last time step when fixed"); }
  arr barV, Vt;
  if(!sys->isKinematic()){
    if(t<T){
      inverse_SymPosDef(barV, Vinv[t+1] + R[t+1]);
      //cout <<"R[t+1]=" <<R[t+1] <<"Vinv[t+1]=" <<Vinv[t+1] <<"barV=" <<barV <<endl;
      Vt = Q[t];
      Vt += B[t]*Hinv[t]*tB[t];
      Vt += barV;
      Vt = Ainv[t]*Vt*invtA[t];
      v[t] = Ainv[t]*(-a[t] + barV*(Vinv[t+1]*v[t+1] + r[t+1]));
      inverse_SymPosDef(Vinv[t](), Vt);
    }
    if(t==T){  //last time slice
      if(!useBwdMsg){
        v[t] = b[t]; //alternative: qhat
        Vinv[t].setDiag(1e-4); //regularization, makes eq (*) above robust
      }else{
        v[T] = bwdMsg_v;
        Vinv[T] = bwdMsg_Vinv;
      }
    }
  }else{
    if(t<T){
      inverse_SymPosDef(barV, Vinv[t+1] + R[t+1]);   //eq (*)
      v[t] = barV * (Vinv[t+1]*v[t+1] + r[t+1]);
      Vt = Hinv[t] + barV;
      inverse_SymPosDef(Vinv[t](), Vt);
    }
    if(t==T){ //last time slice
      if(!useBwdMsg){
        v[t] = b[t]; //alternatives: qhat or b
        Vinv[t].setDiag(1e-0); //regularization, makes eq (*) above robust
      }else{
        v[T] = bwdMsg_v;
        Vinv[T] = bwdMsg_Vinv;
      }
    }
  }
}

void sAICO::updateTaskMessage(uint t, arr& xhat_t){
  
  if(maxStepSize>0. && norm(xhat_t-xhat[t])>maxStepSize){
    arr Delta = xhat_t-xhat[t];
    Delta *= maxStepSize/norm(Delta);
    xhat_t = xhat[t] + Delta;  //really change the given xhat_t (often the belief!!)
  }
  xhat[t]() = xhat_t;
  countSetq++;
  sys->setx(xhat[t]);
  
  //get system matrices
  sys->getControlCosts(NoArr, Hinv[t](), t);
  sys->getDynamics(A[t](), tA[t](), Ainv[t](), invtA[t](), a[t](), B[t](), tB[t](), Q[t](), t);
  sys->getTaskCosts(R[t](), r[t](), t, &rhat(t));
  //double C_alt = scalarProduct(R[t], xhat[t], xhat[t]) - 2.*scalarProduct(r[t], xhat[t]) + rhat(t);
  //cout <<t <<' ' <<C <<' ' <<C_alt <<endl;
}

void sAICO::updateBelief(uint t){
  if(damping && dampingReference.N){
    Binv[t] = Sinv[t] + Vinv[t] + R[t] + damping*eye(R.d1);
    lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t] + damping*dampingReference[t]);
  }else{
    Binv[t] = Sinv[t] + Vinv[t] + R[t];
    lapack_Ainv_b_sym(b[t](), Binv[t], Sinv[t]*s[t] + Vinv[t]*v[t] + r[t]);
  }
}

void sAICO::unrollXhatFwd(uint t){
  CHECK(t>0, "don't update fwd for first time step");
  arr St;
  if(!sys->isKinematic()){
    St = Q[t-1];
    St += B[t-1]*Hinv[t-1]*tB[t-1];
    s[t] = a[t-1] + A[t-1]*xhat[t-1];
    inverse_SymPosDef(Sinv[t](), St);
    
    updateBelief(t);
    xhat[t]() = b[t];
  }else{
    NIY
  }
}