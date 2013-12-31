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
		//cplex.exportModel( "model.lp" );
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
    printX();
		cout << "Solution is " << (isTree() ? "valid" : "invalid") << "\n\n";
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
	IloIntVarArray f(env,2*m);
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "f_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		f[i]=IloIntVar(env,myname.str().c_str());
		
		stringstream myname2;
		myname2 << "f_" << instance.edges.at(i).v2 << "," << instance.edges.at(i).v1;
		f[i+m]=IloIntVar(env,myname2.str().c_str());
	}
	//v
	IloBoolVarArray v(env,n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "v_" << i;
		v[i] = IloBoolVar(env, myname.str().c_str());
	}

	//objective function (1)
	IloExpr objfunc(env);
	for (int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
		objfunc +=x[i+m]*c[i];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

  //(2)
	IloExpr co2(env);
	for(int i=0;i<m;i++){
		co2 += x[i] + x[i+m];
	}
  model.add(co2 == k);
  co2.end();

  //(3)
	IloExpr co3(env);
	for(int i=1;i<n;i++){
		co3 += v[i];
	}
  model.add(co3 == k);
  co3.end();

  //(4)
	IloExpr co4(env); 
	for (int i=0;i<m;i++){ 
		if(instance.edges.at(i).v1==0){ 
			co4+=x[i]; 
			co4+=x[i+m];
		}
	}
	model.add(co4 == 1);
	co4.end();

  //(5) and (6)
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

  //(7)
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

  //(8)
	for (int i=0;i<m;i++){ 
    IloExpr co8_1 = f[i];
    IloExpr co8_2 = f[i+m];
    model.add(0 <= co8_1 <= k);
    model.add(0 <= co8_2 <= k);

    co8_1.end();
    co8_2.end();
	}

  //(9)
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

//  //(10)
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
		model.add(co10_1 == IloMin(1,co10_2));

		co10_1.end();
		co10_2.end();
	}

  // (11)
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

  cout << model << endl;
}

void kMST_ILP::modelMCF()
{
	// ++++++++++++++++++++++++++++++++++++++++++
	// TODO build multi commodity flow model
	// ++++++++++++++++++++++++++++++++++++++++++
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
	IloIntVarArray u(env, n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "u_" << i;
		u[i]=IloIntVar(env, myname.str().c_str());

	}
	//v
	IloBoolVarArray v(env,n);
	for(int i=0;i<n;i++){
		stringstream myname;
		myname << "v_" << i;
		v[i] = IloBoolVar(env, myname.str().c_str());
	}

	//objective function (1)
	IloExpr objfunc(env);
	for (int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
		objfunc +=c[i+m]*x[i+m];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

  //(2)
	for (int i=0;i<m;i++){
		model.add(x[i]+u[instance.edges.at(i).v1] <= u[instance.edges.at(i).v2] + k*(1-x[i]));

		model.add(x[i+m]+u[instance.edges.at(i).v2] <= u[instance.edges.at(i).v1] + k*(1-x[i+m]));
	}
	//(3)
	model.add(u[0]==0);
//(4) 
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
	//(5)	
	for (int i=1;i<n;i++){
		model.add(0<=u[i]);
		model.add(u[i]<=k);
	}

	//(6)
	IloExpr co6(env);
	for (int i=0;i<2*m;i++){
				co6+=x[i];
	}
	model.add(co6 == k);
	co6.end();

	//(7)
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
	//(8)
	IloExpr co8(env);
	for(int i=0;i<n;i++){
		co8+=v[i];
	}
	model.add(co8==k);
	co8.end();
	//(9)
	for(int i=0;i<n;i++){
		IloExpr co9(env);
		co9=u[i];
		model.add(co9<=k*v[i]);
		co9.end();
	}
	//(10)-(11)
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
	//(12)
	IloExpr co12(env);
	for (int i=0;i<n;i++){
		co12+=u[i];
	}
  

  cout << model << endl;
	model.add(co12 == (k*(k+1))/2);
	co12.end();
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
		mat[from][to] = cplex.getValue(x[i]);
		mat[to][from] = cplex.getValue(x[i+m]);
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
