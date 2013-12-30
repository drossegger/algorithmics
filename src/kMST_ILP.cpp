#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance& _instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k ), x(env,instance.n_edges)
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
	// ++++++++++++++++++++++++++++++++++++++++++
	// TODO build single commodity flow model
	// ++++++++++++++++++++++++++++++++++++++++++
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
	for(u_int i=0;i<m;i++){
		stringstream myname;
		myname << "x_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		x[i]=IloBoolVar(env,myname.str().c_str());
	}
	//c
	int c[m];
	for(u_int i=0;i<m;i++){
		c[i]=instance.edges.at(i).weight;
	}
	//u
	IloIntVarArray u(env, n);
	for(u_int i=0;i<n;i++){
		stringstream myname;
		myname << "u_" << i;
		u[i]=IloIntVar(env, myname.str().c_str());

	}

IloExpr myExpr(env);
myExpr += x[1];
myExpr += x[2];
model.add(myExpr >= 1);
myExpr.end(); // IMPORTANT to free memory

	//objective function
	IloExpr objfunc(env);
	for (u_int i=0;i<m; i++){
		objfunc +=x[i]*c[i];
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();

}


bool kMST_ILP::isTreeHelper(int node, bool*& visited, bool**& mat) {
  visited[node] = true;
  ++validEdgeCounter;

  bool ret = true;

  for (u_int i = 0; i < m and ret; i++) {
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
  for(u_int i = 0; i < m; ++i) {
    mat[i] = new bool[m];
  }
  bool* visited = new bool[m];

	for(u_int i=0;i<m;i++){
    int from = instance.edges.at(i).v1;
    int to = instance.edges.at(i).v2;
		mat[from][to] = cplex.getValue(x[i]);
	}

  validEdgeCounter = 0;
  bool valid = isTreeHelper(0, visited, mat); 

  valid = valid and (validEdgeCounter == k+1);

  //clean up
  delete[] visited;
  for(u_int i = 0; i < m; ++i) {
    delete [] mat[i];
  }
  delete [] mat;
   
	return validEdgeCounter;
}

kMST_ILP::~kMST_ILP()
{
	// free global CPLEX resources
	cplex.end();
	model.end();
	env.end();
}
