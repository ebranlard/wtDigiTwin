function B=fCreateBodyRigid(Name,Mass,J_G,rho_G)

B=Body('Name',Name,'Type','Rigid');
B.s_G_inB = rho_G;
B.J_G_inB = J_G  ;
B.Mass    = Mass ;

% B.s_C0_inB = [L; 0; 0];
B.init();
