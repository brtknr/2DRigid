#include <math.h>
#include "Constants.h"

void Invert(int n, double Matrix[][MaxDoF + 1], double Identity[][MaxDoF + 1], double Check[][MaxDoF + 1])
{
    double Reference[MaxDoF + 1][MaxDoF + 1];
    
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++)
            // This is so that Matrix[i][j] is returned unchanged
            Reference[i][j] = Matrix[i][j];	
    }
    
    // Make the leading diagonals 1 and the rest 0 in the inversion matrix
    for (int i = 0; i <= n; i++) {
    	for (int j = 0; j <= n; j++) {
            Identity[i][j] = 0.0;
            Check[i][j] = 0.0;
        }
        Identity[i][i] = 1.0;
    }
    
    for (int i = 0; i <= n; i++) {
    	if (i != n) {
    	    int biggest = i;
    	    for (int k = i + 1; k <= n; k++) {
                if (fabs(Reference[k][i]) > fabs(Reference[biggest][i])) biggest = k;
    	    }
    	    for (int j = 0; j <= n; j++) {
        		double tempdouble;
        		tempdouble = Reference[i][j];
        		Reference[i][j] = Reference[biggest][j];
        		Reference[biggest][j] = tempdouble;
        		tempdouble = Identity[i][j];
        		Identity[i][j] = Identity[biggest][j];
        		Identity[biggest][j] = tempdouble;
 	        }
        }
        
    	double divide = Reference[i][i];
    	for (int j = 0; j <= n; j++) {
    	    Reference[i][j] = Reference[i][j] / divide;
    	    Identity[i][j] = Identity[i][j] / divide;
    	}
        
    	if (i != n) {
    	    for (int k = i + 1; k <= n; k++) {
         		double mult = Reference[k][i];
        		for (int j = 0; j <= n; j++) {
        		    Reference[k][j] -= mult * Reference[i][j];
        		    Identity[k][j] -= mult * Identity[i][j];
        		}
    	    }
    	}
    }
    
    for (int i = n; i >= 1; i -= 1) {
    	for (int k = i - 1; k >= 0; k -= 1) {
    	    double mult = Reference[k][i];
    	    for (int j = 0; j <= n; j++) {
        		Reference[k][j] -= mult * Reference[i][j];
        		Identity[k][j] -= mult * Identity[i][j];
    	    }
    	}
    }
    
    for (int i = 0; i <= n; i++) {
    	for (int j = 0; j <= n; j++) {
        	for (int k = 0; k <= n; k++) {
                Check[i][k] += Matrix[i][j] * Identity[j][k];                	
        	}        	
    	}
    }    
}
