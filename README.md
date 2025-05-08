# initialProfileFoam-v1.0
Custom blastFoam package to set equation-driven initial fields 

Fields are set using equation-defined regions in the createFields.H file. Make changes to this file according to the problem in question.
After editing createFields.H, use wclean and wmake to re-install the package with the new equations defining initial fields.
Example createFields files are located in the History directory.

To apply the fields to a test case, first generate the mesh from the desired test case (e.g. using blockMesh).
Then copy the 0, constant and system files from the desired test case here into initialProfileFoam.
Next, in the 0 directory, rename the field names to match the initialProfileFoam conventions:

p
U
T
alpha1
rho1
rho2 

Where 1 represents "fluid" and 2 represents "gas".
Now run initialProfileFoam in your terminal to apply the fields. You may wish to check the fields by opening open.foam in paraview.
Finally, copy the contents of the 0 directory into the test case and rename them back to the original labelling. You will need to open the files and manually change the 'object' parameter in the OpenFOAM header.

Once this is complete, run your test case as usual in the test directory.
