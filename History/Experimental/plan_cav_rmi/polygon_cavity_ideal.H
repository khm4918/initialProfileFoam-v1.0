 
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

 //low order initialization
          forAll(p,cellI){
    scalar xx=mesh.C()[cellI].x();
    scalar yy=mesh.C()[cellI].y();
	
//fluid boundary geometry params
	scalar y_pos=0;
	scalar lam=0.9e-03;
	scalar a0=0.1e-03;

//shock geometry params
	scalar y_top=-3.9e-03;
	scalar y_bot=-4.1e-03;
	
//cavity size params- number and width
	scalar n = 2;
	scalar w = 1.8e-03;	
	
//polygon vertex positions for cavity zero. nw/2 shift is to find cavity farthest left within boundaries
	std::vector<float> x_p0={0-(n*w)/2};
	scalar y_p0=-0.5e-03;
	
	std::vector<float> x_p1={0.3e-03-(n*w)/2};
	scalar y_p1=-0.7e-03;
	
	std::vector<float> x_p2={0.6e-03-(n*w)/2};
	scalar y_p2=-0.96e-03;
	
	std::vector<float> x_p3={0.9e-03-(n*w)/2};
	scalar y_p3=-1.5e-03;
	
	std::vector<float> x_p4={-0.9e-03-(n*w)/2};
	scalar y_p4=y_p3;
	
	std::vector<float> x_p5={-0.6e-03-(n*w)/2};
	scalar y_p5=y_p2;
	
	std::vector<float> x_p6={-0.3e-03-(n*w)/2};
	scalar y_p6=y_p1;
	
//create vector array to store straight line equations
	std::vector<std::vector<std::vector<float>>> lines; 
	
//calculate subsequent cavitiy positions via translation
	for (int i=0; i<int(n+1); i++) {
		
		x_p0.push_back(x_p0[i]+w);
		x_p1.push_back(x_p1[i]+w);
		x_p2.push_back(x_p2[i]+w);
		x_p3.push_back(x_p3[i]+w);
		x_p4.push_back(x_p4[i]+w);
		x_p5.push_back(x_p5[i]+w);
		x_p6.push_back(x_p6[i]+w);
		
//calculate equations of lines between neighbouring points				
			lines.push_back({
					{(y_p1-y_p0)/(x_p1[i]-x_p0[i]),y_p1-(y_p1-y_p0)/(x_p1[i]-x_p0[i])*x_p1[i]}, //L0
					{(y_p2-y_p1)/(x_p2[i]-x_p1[i]),y_p2-(y_p2-y_p1)/(x_p2[i]-x_p1[i])*x_p2[i]}, //L1
					{(y_p3-y_p2)/(x_p3[i]-x_p2[i]),y_p3-(y_p3-y_p2)/(x_p3[i]-x_p2[i])*x_p3[i]}, //L2
					{(y_p5-y_p4)/(x_p5[i]-x_p4[i]),y_p5-(y_p5-y_p4)/(x_p5[i]-x_p4[i])*x_p5[i]}, //L3
					{(y_p6-y_p5)/(x_p6[i]-x_p5[i]),y_p6-(y_p6-y_p5)/(x_p6[i]-x_p5[i])*x_p6[i]}, //L4
					{(y_p0-y_p6)/(x_p0[i]-x_p6[i]),y_p0-(y_p0-y_p6)/(x_p0[i]-x_p6[i])*x_p0[i]}  //L5
					});
	
	}

					//perturbed fluid-gas boundary and cavities
                    if((yy>=y_pos+a0*(Foam::cos((2*M_PI)/lam*xx)))
					
					//cavity base y position
					|| (yy>=y_p3
					
					and (
					
						//cavity zero
						   (yy<=lines[0][0][0]*xx+lines[0][0][1] and yy<=lines[0][1][0]*xx+lines[0][1][1] and yy<=lines[0][2][0]*xx+lines[0][2][1] and yy<=lines[0][3][0]*xx+lines[0][3][1] and yy<=lines[0][4][0]*xx+lines[0][4][1] and yy<=lines[0][5][0]*xx+lines[0][5][1])
					
						//cavity one
						|| (yy<=lines[1][0][0]*xx+lines[1][0][1] and yy<=lines[1][1][0]*xx+lines[1][1][1] and yy<=lines[1][2][0]*xx+lines[1][2][1] and yy<=lines[1][3][0]*xx+lines[1][3][1] and yy<=lines[1][4][0]*xx+lines[1][4][1] and yy<=lines[1][5][0]*xx+lines[1][5][1])
				
						//cavity two
						|| (yy<=lines[2][0][0]*xx+lines[2][0][1] and yy<=lines[2][1][0]*xx+lines[2][1][1] and yy<=lines[2][2][0]*xx+lines[2][2][1] and yy<=lines[2][3][0]*xx+lines[2][3][1] and yy<=lines[2][4][0]*xx+lines[2][4][1] and yy<=lines[2][5][0]*xx+lines[2][5][1])
						
						)					
				
						
					))
					
					{    //Foam::cos()

                       alpha1[cellI]=0;
					  
					}else{
						
						alpha1[cellI]=1.0;
					
					//shock annulus in water
                    }if((y_bot<=yy)and(yy)<=y_top){    //Foam::cos()

						p[cellI]=8.15e08;
						rho1[cellI]=1200;						
						U[cellI].y()=370.8;

                       }else{
                                           
                       p[cellI]=1.01e5;
					   rho1[cellI]=1000;
                       U[cellI].y()=0.0;
					   
                       }
					   
                      U[cellI].x()=0.0;        
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
    

