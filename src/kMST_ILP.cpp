#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance& _instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), x(env,2*instance.n_edges)
{
	n = instance.n_nodes;
	m = instance.n_edges;
	if( k == 0 ) k = n;
}

void kMST_ILP::solve()
{
	try {
		// initialize CPLEX
		env = IloEnv();
		model = IloModel( env );

		// add model-specific constraints
		if( model_type == "scf" ) modelSCF();
		else if( model_type == "mcf" ) modelMCF();
		else if( model_type == "mtz" ) modelMTZ();
		else {
			cerr << "No existing model chosen\n";
			exit( -1 );
		}

		// build model
		cplex = IloCplex( model );
		// export model to a text file
		cplex.exportModel( "model.lp" );
		// set parameters
		setCPLEXParameters();

		// solve model
		cout << "Calling CPLEX solve ...\n";
		cplex.solve();
		cout << "CPLEX finished.\n\n";
		cout << "CPLEX status: " << cplex.getStatus() << "\n";
		cout << "Branch-and-Bound nodes: " << cplex.getNnodes() << "\n";
		cout << "Objective value: " << cplex.getObjValue() << "\n";
		cout << "CPU time: " << Tools::CPUtime() << "\n\n";
    //printX();
		cout << "Solution is " << (isTree() ? "valid" : "invalid") << "\n\n";
		cplex.writeSolution("solution.lp");
	}
	catch( IloException& e ) {
		cerr << "kMST_ILP: exception " << e << "\n";
		exit( -1 );
	}
	catch( ... ) {
		cerr << "kMST_ILP: unknown exception.\n";
		exit( -1 );
	}
}

// ----- private methods -----------------------------------------------

void kMST_ILP::setCPLEXParameters()
{
	// print every x-th line of node-log and give more details
	cplex.setParam( IloCplex::MIPInterval, 1 );
	cplex.setParam( IloCplex::MIPDisplay, 2 );
	// only use a single thread
	cplex.setParam( IloCplex::Threads, 1 );
}

void kMST_ILP::modelSCF()
{
	//x
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "x_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		x[i]=IloBoolVar(env,myname.str().c_str());
		
		stringstream myname2;
		myname2 << "x_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1;
		x[i+m]=IloBoolVar(env,myname2.str().c_str());
	}
	//c
	int c[2*m];
	for(int i=0;i<m;i++){
		c[i]=instance.edges.at(i).weight;

		c[i+m]=instance.edges.at(i).weight;
	}
	//f
	IloNumVarArray f(env,2*m);
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "f_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		f[i]=IloNumVar(env,myname.str().c_str());
		
		stringstream myname2;
		myname2 << "f_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1;
		f[i+m]=IloNumVar(env,myname2.str().c_str());
	}
	//v
	IloBoolVarArray v(env,n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "v_" << i;
		v[i] = IloBoolVar(env, myname.str().c_str());
	}

	//objective function (1.1)
	IloExpr objfunc(env);
	for (int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
		objfunc +=x[i+m]*c[i];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

  //(3.2)
	IloExpr co2(env);
	for(int i=0;i<m;i++){
		co2 += x[i] + x[i+m];
	}
  model.add(co2 == k);
  co2.end();

  //(3.3)
	IloExpr co3(env);
	for(int i=1;i<n;i++){
		co3 += v[i];
	}
  model.add(co3 == k);
  co3.end();

  //(3.4)
	IloExpr co4(env); 
	for (int i=0;i<m;i++){ 
		if(instance.edges.at(i).v1==0 || instance.edges.at(i).v2==0){ 
			co4+=x[i]; 
			co4+=x[i+m];
		}
	}
	model.add(co4 == 1);
	co4.end();

  //(3.5) and (3.6)
	for (int k=0;k<m;k++){ 
    const int i = instance.edges.at(k).v1;
    const int j = instance.edges.at(k).v2;
    IloExpr co5_1 = x[k];
    IloExpr co5_2 = v[i];
    model.add(co5_1 <= co5_2);
    co5_1.end();
    co5_2.end();

    IloExpr co6_1 = x[k];
    IloExpr co6_2 = v[j];
    model.add(co6_1 <= co6_2);
    co6_1.end();
    co6_2.end();
	}

  //(3.7)
	for (int k=0;k<m;k++){ 
    const int i = instance.edges.at(k).v1;
    const int j = instance.edges.at(k).v2;
    IloExpr co7_1 = v[i];
    co7_1 += x[k];
    co7_1 += x[k+m];

    IloExpr co7_2 = v[j];
    co7_2 += 1;
    model.add(co7_1 <= co7_2);

    co7_1.end();
    co7_2.end();
	}

  //(3.8)
	for (int i=0;i<m;i++){ 
    IloExpr co8_1 = f[i];
    IloExpr co8_2 = f[i+m];
    model.add(0 <= co8_1 <= k);
    model.add(0 <= co8_2 <= k);

    co8_1.end();
    co8_2.end();
	}

  //(3.9)
	IloExpr co9(env); 
	for (int i=0;i<m;i++){ 
		if(instance.edges.at(i).v1==0){ 
			co9+=f[i]; 
		}
		if(instance.edges.at(i).v2==0){
			co9+=f[i+m];
		}
	}
	model.add(co9 == k);
	co9.end();

  //(3.10)
  for(int i=1; i < n; i++) {
		IloNumExpr co10_1(env);
		IloNumExpr co10_2(env);

		for(auto it=instance.incidentEdges.at(i).begin(); it != instance.incidentEdges.at(i).end(); it++) {
			if(instance.edges.at(*it).v1 == i) {	// outgoing edge
				co10_1 -= f[*it];
				co10_1 += f[(*it)+m];

				co10_2 += f[(*it)+m];
			} else {	// incoming edge
				co10_1 += f[*it];
				co10_1 -= f[(*it)+m];

				co10_2 += f[*it];
			}
		}
		model.add(co10_1 == v[i]);

		co10_1.end();
		co10_2.end();
	}

  // (3.11)
  for (int i=0;i<m;i++) { 
    if(instance.edges.at(k).v1 == instance.edges.at(k).v2) {
      continue;
    }

    IloExpr co11_1 = f[i];
    IloExpr co11_2= f[i+m];

    model.add(co11_1 <= k * x[i]);
    model.add(co11_2 <= k * x[i+m]);

    co11_1.end();
    co11_2.end();
  }

	// add
	/*
	for (int i=0;i<m;i++){
		IloExpr coAdd(env);
		coAdd +=x[i];
		coAdd +=x[i+m];

		model.add(coAdd <= 1);
		coAdd.end();

	}
	*/
  //cout << model << endl;
}

void kMST_ILP::modelMCF()
{
	//x
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "x_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		x[i]=IloBoolVar(env,myname.str().c_str());
		
		stringstream myname2;
		myname2 << "x_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1;
		x[i+m]=IloBoolVar(env,myname2.str().c_str());
	}
	//c
	int c[2*m];
	for(int i=0;i<m;i++){
		c[i]=instance.edges.at(i).weight;

		c[i+m]=instance.edges.at(i).weight;
	}
	//f (4.1)
	IloArray<IloBoolVarArray > f(env,n);
	//IloIntVar f[n][2*m];
	for(int l=1;l<n;l++){
		f[l]=IloBoolVarArray(env,2*m);
		for(int i=0;i<m;i++){
			stringstream myname;
			myname << "f_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2<<"^"<<l;
			f[l][i]=IloBoolVar(env,myname.str().c_str());
		
			stringstream myname2;
			myname2 << "f_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1<<"^"<<l;
			f[l][i+m]=IloBoolVar(env,myname2.str().c_str());
		}
	}
	//v
	IloBoolVarArray v(env,n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "v_" << i;
		v[i] = IloBoolVar(env, myname.str().c_str());
	}

	//objective function (1.1)
	IloExpr objfunc(env);
	for (int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
		objfunc +=x[i+m]*c[i];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

  //(4.2)
	IloExpr co2(env);
	for(int i=0;i<m;i++){
		co2 += x[i] + x[i+m];
	}
  model.add(co2 == k);
  co2.end();

  //(4.3)
	IloExpr co3(env);
	for(int i=1;i<n;i++){
		co3 += v[i];
	}
  model.add(co3 == k);
  co3.end();

  //(4.4),(4.5)
	IloExpr co4(env); 
	IloExpr co5(env);
	for (int i=0;i<m;i++){ 
		if(instance.edges.at(i).v1==0){ 
			co4+=x[i]; 
			co5+=x[i+m];
		}
		else if (instance.edges.at(i).v2==0){
			co4+=x[i+m];
			co5+=x[i];
		}
	}
	model.add(co4 == 1);
	model.add(co5==0);
	co4.end();
	co5.end();

	

 //(4.6),(4.7)
	for (int k=0;k<m;k++){ 
    const int i = instance.edges.at(k).v1;
    const int j = instance.edges.at(k).v2;
    IloExpr co5_1 = x[k];
    IloExpr co5_2 = v[i];
    model.add(co5_1 <= co5_2);
    co5_1.end();
    co5_2.end();

    IloExpr co6_1 = x[k];
    IloExpr co6_2 = v[j];
    model.add(co6_1 <= co6_2);
    co6_1.end();
    co6_2.end();
	}
  //(4.8)
	for (int k=0;k<m;k++){ 
    const int i = instance.edges.at(k).v1;
    const int j = instance.edges.at(k).v2;
    IloExpr co7_1 = v[i];
    co7_1 += x[k];
    co7_1 += x[k+m];

    IloExpr co7_2 = v[j];
    co7_2 += 1;
    model.add(co7_1 <= co7_2);

    co7_1.end();
    co7_2.end();
	}
	//(4.9)(4.11)
	for (int l=1;l<n; l++){
		IloExpr co3_9(env);
		IloExpr co3_10(env);
		for (int i=0;i<m;i++){ 
			//(4.9)
			if(instance.edges.at(i).v1==0 ){ 
				co3_9+=f[l][i]; 
			}
			if(instance.edges.at(i).v2==0){
				co3_9+=f[l][i+m];
			}
			//(4.11)
			if(instance.edges.at(i).v2==l && instance.edges.at(i).v1!=l){
				co3_10+=f[l][i];
			}
			else if (instance.edges.at(i).v1==l && instance.edges.at(i).v2!=l){
				co3_10+=f[l][i+m];
			}
    }
		model.add(co3_9 == v[l]);
		//model.add(co3_10==co3_9);
		model.add(co3_10==v[l]);
		co3_9.end();
		co3_10.end();
	}
	//(4.12)
	for (int l=1; l<n;l++){
		for(int i=0;i<m;i++){
			IloExpr co3_11_0(env);
			IloExpr co3_11_3(env);
			IloExpr co3_11_4(env);
			co3_11_3=x[i];
			co3_11_0=f[l][i];
			model.add(0<=co3_11_0);
			model.add(co3_11_0 <=co3_11_3);
			co3_11_0.end();

			co3_11_4=x[i+m];
			IloExpr co3_11_1(env);
			co3_11_1=f[l][i+m];
			model.add(0<=co3_11_1);
			model.add(co3_11_1 <=co3_11_4);
			co3_11_1.end();
			co3_11_3.end();
			co3_11_4.end();
		}
	}
	//(4.13)
	for (int l=1; l<n;l++){
		for (int j=1;j<n;j++){
			IloExpr co33(env);
			bool dontadd=false;
			for(list<u_int>::iterator it=instance.incidentEdges.at(j).begin();it!=instance.incidentEdges.at(j).end();it++){
				if(instance.edges.at(*it).v1==j && j!=l){
					co33+=f[l][*it+m];
					co33-=f[l][*it];
				}
				else if(instance.edges.at(*it).v2==j && j!=l){
					co33+=f[l][*it];
					co33-=f[l][*it +m ];
				}
				else
					dontadd=true;

			}
			if(!dontadd) {
        model.add(co33==0);
      }
			co33.end();

		}
	}

	//(4.14)
	for (int l=1;l<n;l++){
		IloExpr co34(env);
		for (list<u_int>::iterator it=instance.incidentEdges.at(l).begin();it!=instance.incidentEdges.at(l).end();it++){
			if(instance.edges.at(*it).v1==l){
				co34+=f[l][*it];
			}
			else if(instance.edges.at(*it).v2==l){
				co34+=f[l][*it+m];
			}
		}
		model.add(co34 == 0);
		co34.end();
	}


//	//()
//	for (int l=1;l<n;l++){
//		IloExpr co34(env);
//		for (int i=0;i<m;i++){
//			if(instance.edges.at(i).v1==l)
//				co34+=f[l][i+m];
//			if(instance.edges.at(i).v2==l)
//				co34+=f[l][i];
//		}
//		model.add(co34-v[l]==0);
//		co34.end();
//	}
	//35
	/*
	for( int l=0;l<n-1;l++){
		IloExpr co35_0(env);
		IloExpr co35_1(env);
		for(int i=0;i<m;i++){
			if(instance.edges.at(i).v2==l && instance.edges.at(i).v1!=l){
				co35_0+=f[l][i];
			}
			else if (instance.edges.at(i).v1==l && instance.edges.at(i).v2!=l){
				co35_0+=f[l][i+m];
			}
		}
		for (int i=0;i<m;i++){ 
			if(instance.edges.at(i).v1==0 ){ 
				co35_1+=f[l][i]; 
			}
			if(instance.edges.at(i).v2==0){
				co35_1+=f[l][i+m];
			}	
		}
		model.add(co35_0 == co35_1);
		co35_0.end();
		co35_1.end();
	}*/
	//(4.16)
	IloExpr co(env);	
	for (int l=1;l<n;l++){
		for(list<u_int>::iterator it=instance.incidentEdges.at(l).begin();it!=instance.incidentEdges.at(l).end();it++){
			if(instance.edges.at(*it).v1==l && instance.edges.at(*it).v2!=l){
				co+=f[l][*it+m];
			}
			else if(instance.edges.at(*it).v2==l && instance.edges.at(*it).v1!=l){
				co+=f[l][*it];
			}
		}
	}
  model.add(co==k);
  co.end();
    
	//(4.10)
	IloExpr co_35(env);
	for (int l=1;l<n;l++){
		for(int j=0;j<m;j++){
			const int v1=instance.edges.at(j).v1;
			const int v2=instance.edges.at(j).v2;
			if(v1 == 0){
				co_35 +=f[l][j];
			}
			else if(v2 == 0){
				co_35 +=f[l][j+m];
			}
		}
	}
	model.add(co_35 == k);
	co_35.end();
  
		
	//(4.15)
	for (int l=1;l<n;l++){
		IloExpr co(env);
		for (list<u_int>::iterator it=instance.incidentEdges.at(0).begin();it!=instance.incidentEdges.at(0).end();it++){
			co+=f[l][*it];
		}
		model.add(co - v[l]==0);
		co.end();
	}
	//(4.17)
  for(int k=0;k<m;k++){
		IloExpr co_0(env);
		IloExpr co_1(env);
    for(int l=1;l<n;l++){
      co_0+=f[l][k];
      co_1+=f[l][k+m];
    }

    IloExpr x_ij(env);
    IloExpr x_ji(env);

    x_ij += x[k];
    x_ji += x[k+m];
		model.add(co_0>=x_ij);
		model.add(co_1>=x_ji);
		co_0.end();
		co_1.end();

    x_ij.end();
    x_ji.end();
	}

  //cout << model;
}

void kMST_ILP::modelMTZ()
{
	//x
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "x_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		x[i]=IloBoolVar(env,myname.str().c_str());
		
		stringstream myname2;
		myname2 << "x_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1;
		x[i+m]=IloBoolVar(env,myname2.str().c_str());
		
	}
	//c
	int c[2*m];
	for(int i=0;i<m;i++){
		c[i]=instance.edges.at(i).weight;

		c[i+m]=instance.edges.at(i).weight;
	}
	//u
	IloNumVarArray u(env, n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "u_" << i;
		u[i]=IloNumVar(env, myname.str().c_str());

	}
	//v
	IloBoolVarArray v(env,n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "v_" << i;
		v[i] = IloBoolVar(env, myname.str().c_str());
	}

	//objective function (1.1)
	IloExpr objfunc(env);
	for (int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
		objfunc +=c[i+m]*x[i+m];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

  //(2.1)
	for (int i=0;i<m;i++){
		model.add(x[i]+u[instance.edges.at(i).v1] <= u[instance.edges.at(i).v2] + k*(1-x[i]));

		model.add(x[i+m]+u[instance.edges.at(i).v2] <= u[instance.edges.at(i).v1] + k*(1-x[i+m]));
	}
	//(2.2)
	model.add(u[0]==0);
  //(2.3) 
	IloExpr co4(env); 
	for (int i=0;i<m;i++){ 
		if(instance.edges.at(i).v1==0 /*or instance.edges.at(i).v2==0*/){ 
			co4+=x[i]; 
		}
		if(instance.edges.at(i).v2==0){
			co4+=x[i+m];
		}
	}
	model.add(co4 == 1);
	co4.end();
	//(2.4)	
	for (int i=1;i<n;i++){
		model.add(0<=u[i]);
		model.add(u[i]<=k);
	}

	//(2.5)
	IloExpr co6(env);
	for (int i=0;i<2*m;i++){
				co6+=x[i];
	}
	model.add(co6 == k);
	co6.end();

	//(2.6)
	for(int i=1;i<n;i++){
		IloExpr co7(env);
		for(list<u_int>::iterator it=instance.incidentEdges.at(i).begin();it != instance.incidentEdges.at(i).end();++it){
			if(instance.edges.at(*it).v2==i)
				co7+=x[*it];
			else
				co7+=x[*it+m];
		}
		model.add(co7 <= 1);
		co7.end();
	}
	//(2.7)
	IloExpr co8(env);
	for(int i=0;i<n;i++){
		co8+=v[i];
	}
	model.add(co8==k);
	co8.end();
	//(2.8)
	for(int i=0;i<n;i++){
		IloExpr co9(env);
		co9=u[i];
		model.add(co9<=k*v[i]);
		co9.end();
	}
	//(2.9)-(2.10)
	for(int i=0;i<m;i++){
			if(instance.edges.at(i).v1!=0 and instance.edges.at(i).v2!=0){
				IloExpr co10(env),co11(env),co11_1(env),co10_1(env);
				co10=x[i];
				co11=x[m+i];
				co10_1=u[instance.edges.at(i).v1];
				co11_1=u[instance.edges.at(i).v2];
				model.add(co10<=co10_1);
				model.add(co11<=co10_1);
				model.add(co10<=co11_1);
				model.add(co11<=co11_1);
				co10.end();co11.end();co11_1.end();co10_1.end();
			}
	}
	//(2.11)
	IloExpr co12(env);
	for (int i=0;i<n;i++){
		co12+=u[i];
	}
	model.add(co12 == (k*(k+1))/2);
	co12.end();

	//2.12
	for (int i=0;i<m;i++){
	IloExpr coAdd(env);
	coAdd+=x[i];
	coAdd+=x[i+m];
	model.add(coAdd <=1);
	coAdd.end();
	}
}


bool kMST_ILP::isTreeHelper(int node, bool*& visited, bool**& mat) {
  visited[node] = true;
  ++validEdgeCounter;

  bool ret = true;

  for (int i = 0; i < m and ret; i++) {
    if(i == node) continue;

    if(mat[node][i]) {
      if(visited[i]) {
        return false;
      }
      ret = ret and isTreeHelper(i, visited, mat);
    }
  }
   
  return ret;
}

bool kMST_ILP::isTree() {
  //helping variables
  bool** mat = new bool*[m];
  for(int i = 0; i < m; ++i) {
    mat[i] = new bool[m];
    memset(mat[i], 0, sizeof(bool) * m);
  }

  bool* visited = new bool[m];
  memset(visited, 0, sizeof(bool) * m);

	for(int i=0;i<m;i++){
    int from = instance.edges.at(i).v1;
    int to = instance.edges.at(i).v2;
		mat[from][to] = 0.9 < cplex.getValue(x[i]);
		mat[to][from] = 0.9 < cplex.getValue(x[i+m]);
	}

  validEdgeCounter = 0;
  bool valid = isTreeHelper(0, visited, mat); 

  valid = valid and (validEdgeCounter == k+1);

  //clean up
  delete[] visited;
  for(int i = 0; i < m; ++i) {
    delete [] mat[i];
  }
  delete [] mat;
   
	return valid;
}

void kMST_ILP::printX() {
  cout << "Decision Variables:" << endl;
	for(int i=0;i<m;i++){
    int from = instance.edges.at(i).v1;
    int to = instance.edges.at(i).v2;
    try{
      cout << "x_" << from << "," << to << " (" << i << "): \t" << cplex.getIntValue(x[i]) << endl;
    }catch (IloAlgorithm::NotExtractedException& e) {
      cout << "x_" << from << "," << to << " (" << i << "): \t not used by solution" << endl;
    }
    try{
      cout << "x_" << to << "," << from << " (" << i << "): \t" << cplex.getIntValue(x[i+m]) << endl;
    }catch (IloAlgorithm::NotExtractedException& e) {
      cout << "x_" << to << "," << from << " (" << i << "): \t not used by solution" << endl;
    }
  }
}

kMST_ILP::~kMST_ILP()
{
  // free global CPLEX resources
  cplex.end();
  model.end();
  env.end();
}
