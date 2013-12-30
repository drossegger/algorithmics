#include "kMST_ILP.h"

kMST_ILP::kMST_ILP( Instance& _instance, string _model_type, int _k ) :
	instance( _instance ), model_type( _model_type ), k( _k )
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
	IloBoolVarArray x(env,m);
	for(int i=0;i<m;i++){
		stringstream myname;
		myname << "x_" << instance.edges.at(i).v1 << "," <<instance.edges.at(i).v2;
		x[i]=IloBoolVar(env,myname.str().c_str());
	}
	//c
	int c[m];
	for(int i=0;i<m;i++){
		c[i]=instance.edges.at(i).weight;
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
	}
	model.add(IloMinimize(env,objfunc));
	objfunc.end();
	
	//(2)
	for (int i=0;i<m;i++){
		IloExpr co3_0(env);
		co3_0+=k*x[i];
		co3_0+=u[instance.edges.at(i).v1];

		IloExpr co3_1(env);
		co3_1+=u[instance.edges.at(i).v2];
		co3_1+=k-1;
		model.add(co3_0 <= co3_1);
		co3_0.end();
		co3_1.end();
	}
	//(3)
	IloExpr co3(env);
	co3=u[0];
	model.add(co3==0);
	co3.end();

	//(4)
	IloExpr co4(env);
	for (int i=0;i<m;i++){
		if(instance.edges.at(i).v1==0 or instance.edges.at(i).v2==0){
			co4+=x[i];
		}
	}
	model.add(co4 == 1);
	co4.end();
	//(5)	
	for (int i=0;i<n;i++){
		model.add(0<=u[i]);
		model.add(u[i]<=k);
	}

	//(6)
	IloExpr co6(env);
	for (int i=0;i<m;i++){
				co6+=x[i];
	}
	model.add(co6 == k);
	co6.end();

	//(7)
	for(int i=0;i<n;i++){
		IloExpr co7(env);
		for(int j=0;j<m;j++){
			if(instance.edges.at(j).v2==i)
				co7+=x[j];
		}
		model.add(co7 <= 1);
		co7.end();
	}
	//(8)
	IloExpr co8(env);
	for(int i=1;i<n;i++){
		co8+=v[i];
	}
	model.add(co8==k);
	co8.end();
	//(9)
	for(int i=1;i<n;i++){
		model.add(u[i]<=k*v[i]);
	}
	//(10)-(11)
	for(int i=0;i<m;i++){
		if(instance.edges.at(i).v1!=0)
		model.add(x[i]<=u[instance.edges.at(i).v1]);
		if(instance.edges.at(i).v2!=0)
		model.add(x[i]<=u[instance.edges.at(i).v2]);
	}
	//(12)
	IloExpr co12(env);
	for (int i=0;i<n;i++){
		co12+=u[i];
	}
	model.add(co12 == (k*(k+1))/2);
	co12.end();
}

kMST_ILP::~kMST_ILP()
{
	// free global CPLEX resources
	cplex.end();
	model.end();
	env.end();
}
