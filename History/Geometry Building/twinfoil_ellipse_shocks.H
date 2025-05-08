 
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
	
//twinfoil one params
	scalar x_twin1=-3.0e-03;
	scalar y_twin1=3.75e-03;
	scalar a_twin1=0.05e-03;
	scalar b_twin1=1.25e-03;
	
//twinfoil two params
	scalar x_twin2=-3.0e-03;
	scalar y_twin2=-3.75e-03;
	scalar a_twin2=0.05e-03;
	scalar b_twin2=1.25e-03;
	
//back foil two params
	scalar x_back=-8.0e-03;
	scalar y_back=0e-03;
	scalar a_back=0.05e-03;
	scalar b_back=2.5e-03;
	
					//twinfoil one
                    if((((xx-x_twin1)*(xx-x_twin1))/(a_twin1*a_twin1)+((yy-y_twin1)*(yy-y_twin1))/(b_twin1*b_twin1))<=1){   //Foam::cos()
					
						if((xx==x_twin1)and(yy==y_twin1)){ 
							p[cellI]=8.15e08;
							rho1[cellI]=1200;						
							U[cellI].x()=0.0;
							U[cellI].y()=0.0;
							
							Info<< "Found discontinuity" << endl;
							}
						
						else{						
							p[cellI]=8.15e08;
							rho1[cellI]=1200;
							U[cellI].x()=370.8/Foam::sqrt(4*((xx-x_twin1)*(xx-x_twin1)/(a_twin1*a_twin1*a_twin1*a_twin1)+(yy-y_twin1)*(yy-y_twin1)/(b_twin1*b_twin1*b_twin1*b_twin1)))*2*(xx-x_twin1)/(a_twin1*a_twin1);
							U[cellI].y()=370.8/Foam::sqrt(4*((xx-x_twin1)*(xx-x_twin1)/(a_twin1*a_twin1*a_twin1*a_twin1)+(yy-y_twin1)*(yy-y_twin1)/(b_twin1*b_twin1*b_twin1*b_twin1)))*2*(yy-y_twin1)/(b_twin1*b_twin1);
							
							}
							
					//twinfoil two						
                    }else if((((xx-x_twin2)*(xx-x_twin2))/(a_twin2*a_twin2)+((yy-y_twin2)*(yy-y_twin2))/(b_twin2*b_twin2))<=1){   //Foam::cos()
					
						if((xx==x_twin2)and(yy==y_twin2)){ 
							p[cellI]=8.15e08;
							rho1[cellI]=1200;						
							U[cellI].x()=0.0;
							U[cellI].y()=0.0;
							
							Info<< "Found discontinuity" << endl;
							}
						
						else{						
							p[cellI]=8.15e08;
							rho1[cellI]=1200;
							U[cellI].x()=370.8/Foam::sqrt(4*((xx-x_twin2)*(xx-x_twin2)/(a_twin2*a_twin2*a_twin2*a_twin2)+(yy-y_twin2)*(yy-y_twin2)/(b_twin2*b_twin2*b_twin2*b_twin2)))*2*(xx-x_twin2)/(a_twin2*a_twin2);
							U[cellI].y()=370.8/Foam::sqrt(4*((xx-x_twin2)*(xx-x_twin2)/(a_twin2*a_twin2*a_twin2*a_twin2)+(yy-y_twin2)*(yy-y_twin2)/(b_twin2*b_twin2*b_twin2*b_twin2)))*2*(yy-y_twin2)/(b_twin2*b_twin2);
							
							}
							
					//back foil					
                    }else if((((xx-x_back)*(xx-x_back))/(a_back*a_back)+((yy-y_back)*(yy-y_back))/(b_back*b_back))<=1){   //Foam::cos()
					
						if((xx==x_back)and(yy==y_back)){ 
							p[cellI]=8.15e08;
							rho1[cellI]=1200;						
							U[cellI].x()=0.0;
							U[cellI].y()=0.0;
							
							Info<< "Found discontinuity" << endl;
							}
						
						else{						
							p[cellI]=8.15e08;
							rho1[cellI]=1200;
							U[cellI].x()=370.8/Foam::sqrt(4*((xx-x_back)*(xx-x_back)/(a_back*a_back*a_back*a_back)+(yy-y_back)*(yy-y_back)/(b_back*b_back*b_back*b_back)))*2*(xx-x_back)/(a_back*a_back);
							U[cellI].y()=370.8/Foam::sqrt(4*((xx-x_back)*(xx-x_back)/(a_back*a_back*a_back*a_back)+(yy-y_back)*(yy-y_back)/(b_back*b_back*b_back*b_back)))*2*(yy-y_back)/(b_back*b_back);
							
							}
							
					   }else{
                                           
                       p[cellI]=1.01e5;
					   rho1[cellI]=1000;
                       U[cellI].x()=0.0;
                       U[cellI].y()=0.0;
					  
                       }
        

                      T[cellI]=298.0;  
                      rho2[cellI]=1.293;
					  alpha1[cellI]=1.0;
       
                     }
      // high order initialization    
      
        /* label np=40;
                    
         forAll(U,cellI)
        {
            scalar sumy_twin1=0.0;
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
                            sumy_twin1+=1.0e-6;
                            sumY1+=1.0e-6;
                            sumY2+=0.999998;
                        }
                        else
                        {
                            sumy_twin1+=0.999998;
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

             Y[0][cellI] = sumy_twin1/sumNum;
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
    

