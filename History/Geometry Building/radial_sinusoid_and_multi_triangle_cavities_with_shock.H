 
 /*   dimensionedScalar R
    (
        transportProperties.lookup("R")
    );

    dimensionedScalar Cv
    (
        transportProperties.lookup("Cv")
    );*/

    Info<< "Reading field pressure\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    

    Info<< "Reading field rho\n" << endl;
    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
        volScalarField rho2
    (
        IOobject
        (
            "rho2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField alpha1    
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field T\n" << endl;
  volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

 Info<< "Finish initialization field" << endl;

//fluid boundary geometry params
	scalar r0=2.5e-03;
	scalar a0=0.1e-03;
	int n = 16;

//shock geometry params
	scalar r_in=4.9e-03;
	scalar r_out=5.1e-03;
	
//first cavity triangle positions
	std::vector<float> x_p0={0};
	std::vector<float> y_p0={2.79e-03};
	
	std::vector<float> x_p1={1.42e-03};
	std::vector<float> y_p1={3.51e-03};
	
	std::vector<float> x_p2={-1.42e-3};
	std::vector<float> y_p2={3.51e-03};
				
	std::vector<std::vector<std::vector<float>>> lines; 

//rotation matrix params
	scalar rot=M_PI*2/(n/2);
				
//subsequent cavity triangle positions

	for (int i=1; i<int(n/2 + 1); i++) {
		
		x_p0.push_back(x_p0[0]*Foam::cos(rot*i)-y_p0[0]*Foam::sin(rot*i));
		y_p0.push_back(x_p0[0]*Foam::sin(rot*i)+y_p0[0]*Foam::cos(rot*i));
		
		x_p1.push_back(x_p1[0]*Foam::cos(rot*i)-y_p1[0]*Foam::sin(rot*i));
		y_p1.push_back(x_p1[0]*Foam::sin(rot*i)+y_p1[0]*Foam::cos(rot*i));
		
		x_p2.push_back(x_p2[0]*Foam::cos(rot*i)-y_p2[0]*Foam::sin(rot*i));
		y_p2.push_back(x_p2[0]*Foam::sin(rot*i)+y_p2[0]*Foam::cos(rot*i));
		
		if(x_p2[i-1]-x_p1[i-1] == 0){
			
			cout << "found floating point exception \n";
			
			lines.push_back({
					{(y_p1[i-1]-y_p0[i-1])/(x_p1[i-1]-x_p0[i-1]),y_p1[i-1]-(y_p1[i-1]-y_p0[i-1])/(x_p1[i-1]-x_p0[i-1])*x_p1[i-1]},
					{0,0},
					{(y_p2[i-1]-y_p0[i-1])/(x_p2[i-1]-x_p0[i-1]),y_p2[i-1]-(y_p2[i-1]-y_p0[i-1])/(x_p2[i-1]-x_p0[i-1])*x_p2[i-1]}
					});	
					
		}
		
		else {
			
			lines.push_back({
					{(y_p1[i-1]-y_p0[i-1])/(x_p1[i-1]-x_p0[i-1]),y_p1[i-1]-(y_p1[i-1]-y_p0[i-1])/(x_p1[i-1]-x_p0[i-1])*x_p1[i-1]},
					{(y_p2[i-1]-y_p1[i-1])/(x_p2[i-1]-x_p1[i-1]),y_p2[i-1]-(y_p2[i-1]-y_p1[i-1])/(x_p2[i-1]-x_p1[i-1])*x_p2[i-1]},
					{(y_p2[i-1]-y_p0[i-1])/(x_p2[i-1]-x_p0[i-1]),y_p2[i-1]-(y_p2[i-1]-y_p0[i-1])/(x_p2[i-1]-x_p0[i-1])*x_p2[i-1]}
					});
					
			}
		}

 //low order initialization
          forAll(p,cellI){
    scalar xx=mesh.C()[cellI].x();
    scalar yy=mesh.C()[cellI].y();					
	
					//radial perturbation to fluid-gas boundary
                    if((xx*xx+yy*yy)<=(r0-a0*Foam::cos(n*Foam::atan(yy/xx)))*(r0-a0*Foam::cos(n*Foam::atan(yy/xx)))
					
					//cavity one
					|| (yy<=lines[0][1][0]*xx+lines[0][1][1] and yy>=lines[0][0][0]*xx+lines[0][0][1] and yy>=lines[0][2][0]*xx+lines[0][2][1])
					
					//cavity two
					|| (yy<=lines[1][1][0]*xx+lines[1][1][1] and yy>=lines[1][0][0]*xx+lines[1][0][1] and yy>=lines[1][2][0]*xx+lines[1][2][1])
					
					//cavity three *vertical line
					|| (yy<=lines[2][0][0]*xx+lines[2][0][1] and xx>=x_p1[2] and yy>=lines[2][2][0]*xx+lines[2][2][1])
				
					//cavity four
					|| (yy<=lines[3][0][0]*xx+lines[3][0][1] and yy>=lines[3][1][0]*xx+lines[3][1][1] and yy<=lines[3][2][0]*xx+lines[3][2][1])
					
					//cavity five	
					|| (yy<=lines[4][0][0]*xx+lines[4][0][1] and yy>=lines[4][1][0]*xx+lines[4][1][1] and yy<=lines[4][2][0]*xx+lines[4][2][1])
					
					//cavity six
					|| (yy<=lines[5][0][0]*xx+lines[5][0][1] and yy>=lines[5][1][0]*xx+lines[5][1][1] and yy<=lines[5][2][0]*xx+lines[5][2][1])
					
					//cavity seven *vertical line
					|| (yy>=lines[6][0][0]*xx+lines[6][0][1] and xx<=x_p1[6] and yy<=lines[6][2][0]*xx+lines[6][2][1])
					
					//cavity eight
					|| (yy>=lines[7][0][0]*xx+lines[7][0][1] and yy<=lines[7][1][0]*xx+lines[7][1][1] and yy>=lines[7][2][0]*xx+lines[7][2][1])
					
					)
						
					{    //Foam::cos()

                       alpha1[cellI]=0;
					  
					}else{
						
						alpha1[cellI]=1.0; 
					
					//shock annulus in water
                    }if((r_in*r_in<=(xx*xx+yy*yy))and(xx*xx+yy*yy)<=r_out*r_out){    //Foam::cos()

						p[cellI]=8.15e08;
						rho1[cellI]=1200;						
						U[cellI].x()=-(370.8/Foam::sqrt(xx*xx+yy*yy))*xx;
						U[cellI].y()=-(370.8/Foam::sqrt(xx*xx+yy*yy))*yy;

                       }else{
                                           
                       p[cellI]=1.01e5;
					   rho1[cellI]=1000;
                       U[cellI].x()=0.0;
                       U[cellI].y()=0.0;
					   
                       }
        
                      T[cellI]=298.0;        
                      rho2[cellI]=1.293;
       
                     }
      // high order initialization    
      
        /* label np=40;
                    
         forAll(U,cellI)
        {
            scalar sumY0=0.0;
            scalar sumY1=0.0;
            scalar sumY2=0.0;
            //vector sumV=vector::zero;
            scalar sumNum=0.0;
            if(mixMeshType[cellI]==0)
            {
                label pt[4];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];
                pt[3]=mixVertLabels[cellI][3];

                scalar xM[4],yM[4];
                for(label m=0;m<4;m++)
                {
                    xM[m]=mesh.points()[pt[m]].x(); 
                    yM[m]=mesh.points()[pt[m]].y();
                }

                scalar N[4];
                for(int i=0; i<=np; i++)
                {
                    for(int j=0; j<=np; j++)
                    {
                        scalar wgI = -1.0+double(i)*2.0/double(np);
                        scalar wgJ = -1.0+double(j)*2.0/double(np);
                        N[0]=(1.0-wgI)*(1.0-wgJ)/4.0;
                        N[1]=(1.0+wgI)*(1.0-wgJ)/4.0;
                        N[2]=(1.0+wgI)*(1.0+wgJ)/4.0;
                        N[3]=(1.0-wgI)*(1.0+wgJ)/4.0;    
                        scalar xx = N[0]*xM[0]+N[1]*xM[1]+N[2]*xM[2]+N[3]*xM[3];
                        scalar yy = N[0]*yM[0]+N[1]*yM[1]+N[2]*yM[2]+N[3]*yM[3];
                        if(((xx-0.5)*(xx-0.5)+(yy-0.5)*(yy-0.5))<=0.0625)
                        {
                            sumY0+=1.0e-6;
                            sumY1+=1.0e-6;
                            sumY2+=0.999998;
                        }
                        else
                        {
                            sumY0+=0.999998;
                            sumY1+=1.0e-6;
                            sumY2+=1.0e-6;
                        }
                        sumNum+=1.0;
                    }
                }
            }
            else
            {
                label pt[3];
                pt[0]=mixVertLabels[cellI][0];
                pt[1]=mixVertLabels[cellI][1];
                pt[2]=mixVertLabels[cellI][2];

                scalar xM[3],yM[3];
                for(label m=0;m<3;m++)
                {
                   xM[m]=mesh.points()[pt[m]].x(); 
                   yM[m]=mesh.points()[pt[m]].y();
                } 

                scalar N[3];

           }

             Y[0][cellI] = sumY0/sumNum;
             Y[1][cellI] = sumY1/sumNum;
             Y[2][cellI] = sumY2/sumNum;
            p[cellI]=1.01e5;
            T[cellI]=298.0;  
           // rho[cellI] = sumR/sumNum;
        }
        */



    U.write();
    p.write();
    alpha1.write();
    rho1.write();
    rho2.write();
    T.write();
    U.write();
    

