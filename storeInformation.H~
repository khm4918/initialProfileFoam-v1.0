  label noe=mesh.nCells();
    labelListList mixFaceLabels(noe);
    labelListList mixVertLabels(noe);
    labelList mixMeshType(noe);

    hexMatcher hex;
    prismMatcher prism;
    cellShape shape;
    const pointField& faceCentres = mesh.faceCentres();
    forAll(mixFaceLabels,cellI)
    {
       bool isHex = hex.matches(mesh,cellI, shape);
       bool isPrism = prism.matches(mesh,cellI, shape);
       if(isHex)
       {
          mixMeshType[cellI]=0;
          label faceNum=hex.faceLabels().size();
          label vertNum=hex.vertLabels().size();
          mixFaceLabels[cellI].setSize(faceNum);
          mixVertLabels[cellI].setSize(vertNum);
          forAll(hex.faceLabels(),cfI) mixFaceLabels[cellI][cfI]=hex.faceLabels()[cfI];
          forAll(hex.vertLabels(),cpI) mixVertLabels[cellI][cpI]=hex.vertLabels()[cpI];

          label pt[8],fe[6];
          forAll(mixVertLabels[cellI],pI) pt[pI]=mixVertLabels[cellI][pI];
          forAll(mixFaceLabels[cellI],fI) fe[fI]=mixFaceLabels[cellI][fI];

          label Face[2]={4,5};
          label i=0;
          for(label m=0; m<6; m++)
          {
             if(mag(faceCentres[fe[m]].z()-mesh.C()[cellI].z())>SMALL)
             {
                if(i>1) 
                {
                   Info<<"wrong face arrangement or not in X-Y plane!!"<<endl;
                   exit(0);
                }
                Face[i++]=m;
             }
          }
    
          if(Face[0]==0 && Face[1]==1)
          {
             mixVertLabels[cellI][0]=pt[0];
             mixVertLabels[cellI][1]=pt[3];
             mixVertLabels[cellI][2]=pt[7];
             mixVertLabels[cellI][3]=pt[4];
             mixVertLabels[cellI][4]=pt[1];
             mixVertLabels[cellI][5]=pt[2];
             mixVertLabels[cellI][6]=pt[6];
             mixVertLabels[cellI][7]=pt[5];
             mixFaceLabels[cellI][0]=fe[2];
             mixFaceLabels[cellI][1]=fe[3];
             mixFaceLabels[cellI][2]=fe[4];
             mixFaceLabels[cellI][3]=fe[5];
             mixFaceLabels[cellI][4]=fe[0];
             mixFaceLabels[cellI][5]=fe[1];
          }
          else if(Face[0]==2 && Face[1]==3)
          {
             mixVertLabels[cellI][0]=pt[0];
             mixVertLabels[cellI][1]=pt[1];
             mixVertLabels[cellI][2]=pt[5];
             mixVertLabels[cellI][3]=pt[4];
             mixVertLabels[cellI][4]=pt[3];
             mixVertLabels[cellI][5]=pt[2];
             mixVertLabels[cellI][6]=pt[6];
             mixVertLabels[cellI][7]=pt[7];
             mixFaceLabels[cellI][0]=fe[0];
             mixFaceLabels[cellI][1]=fe[1];
             mixFaceLabels[cellI][2]=fe[4];
             mixFaceLabels[cellI][3]=fe[5];
             mixFaceLabels[cellI][4]=fe[2];
             mixFaceLabels[cellI][5]=fe[3];
          }
       }
       else if(isPrism)
       {
          mixMeshType[cellI]=1;
          label faceNum=prism.faceLabels().size();
          label vertNum=prism.vertLabels().size();
          mixFaceLabels[cellI].setSize(faceNum);
          mixVertLabels[cellI].setSize(vertNum);
          forAll(prism.faceLabels(),cfI) mixFaceLabels[cellI][cfI]=prism.faceLabels()[cfI];
          forAll(prism.vertLabels(),cpI) mixVertLabels[cellI][cpI]=prism.vertLabels()[cpI];

          label pt[6],fe[5];
          forAll(mixVertLabels[cellI],pI) pt[pI]=mixVertLabels[cellI][pI];
          forAll(mixFaceLabels[cellI],fI) fe[fI]=mixFaceLabels[cellI][fI];

          mixVertLabels[cellI][0]=pt[2];
          mixVertLabels[cellI][1]=pt[0];
          mixVertLabels[cellI][2]=pt[1];
          mixVertLabels[cellI][3]=pt[5];
          mixVertLabels[cellI][4]=pt[3];
          mixVertLabels[cellI][5]=pt[4];
          mixFaceLabels[cellI][0]=fe[2];
          mixFaceLabels[cellI][1]=fe[3];
          mixFaceLabels[cellI][2]=fe[4];
          mixFaceLabels[cellI][3]=fe[0];
          mixFaceLabels[cellI][4]=fe[1];
       }
       else
       {
          Info<<"this scheme only support Tri and Qua mesh element for 2d!!"<<endl;
          exit(0);
       }
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neigh = mesh.neighbour();

    //Boundary Point_Search





