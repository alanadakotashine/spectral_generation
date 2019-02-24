#include "mex.h"
#include <unistd.h>
#include <sys/mman.h>
#include <assert.h>
#include <strings.h>
#include <math.h>
#include "uthash.h"
#include "heap.h"
#include <stdlib.h>   /* malloc */
#include <stdio.h> 




///Structure for hashtable entry. Hashtable contains a struct for each value that can be added to the eigenScore.
//each value can corresond to multiple edges. edgeLocs contains the pairs for corresponding to all edges with score
//equal to value. current stores the index of the edge we are considering. 
typedef struct value_edgeLocs_entry {
    double id;
    int* edgeLocs;
    int size;
    int capacity;
    int current;
    UT_hash_handle hh;
} value_edgeLocs_entry;

// Structure for a single heap entry. Heap is just used for sorting. 
typedef struct heap_entry {
    void* key;   // Key for this entry
    void* value; // Value for this entry
} heap_entry;

//The compare function must take floats, the default compare function in the heap implementation linked in the read me
//compares ints. Note table is a pointer, not a pointer of pointers as in the heap implementation in the read me. 
//The heap implementation in the read me must be changed to store the table as a pointer
typedef struct heap {
    int (*compare_func)(void*, void*); // The key comparison function to use
    int active_entries;  // The number of entries in the heap
    int minimum_pages;   // The minimum number of pages to maintain, based on the initial cap.
    int allocated_pages; // The number of pages in memory that are allocated
    heap_entry* table; 
} heap;

void record(int *data, int i, int j, int k, int l, int entry, int point){
    *(data+(point*12)+(4*entry)+0)=i;
    *(data+(point*12)+(4*entry)+1)=j;
    *(data+(point*12)+(4*entry)+2)=k;
    *(data+(point*12)+(4*entry)+3)=l;


}

double compAlpha(double add1, double add2, double remove1, double remove2){
    double alpha = min(remove1,remove2);
    double alpha2 = min(1-add1,1-add2);
    alpha = min(alpha,alpha2);
    return alpha;
}

void compValueSingleEdge(double *value, double *u, double *v, double *degree, int i, int j, double *A, mwSize n, int AddFlag, double eigNorm, double targetSum){
    double edgeChange;
    if (AddFlag == 1){
        //if i'm addding an edge, I'm adding 1-A_ij. But then it's value is multiplied by -1
        edgeChange = A[i*n+j]-1;
    }
    else{
        //when i remove an edge, should be on opposite sides and i'm removing edge value
        edgeChange = A[i*n+j];
    }

    double normScore = (edgeChange*edgeChange)/(degree[i]*degree[j]);
    //negative in front of additions. this shoud stay negative when i add items on same side of cut
    //double eigScoreAdd = ((-1*alpha*u[i]*v[j])/(sqrt(degree[i])*sqrt(degree[j])))+((-1*alpha*u[k]*v[l])/(sqrt(degree[k])*sqrt(degree[l])));
    //double eigScoreNeg = ((alpha*u[i]*v[k])/(sqrt(degree[i])*sqrt(degree[k])))+((alpha*u[j]*v[l])/(sqrt(degree[j])*sqrt(degree[l])));
    double eigScore = edgeChange*u[i]*v[j];
    //mexPrintf("eigScore %f\n", eigScore);
    *(value) = eigScore/(eigNorm*sqrt(degree[i])*sqrt(degree[j]));
    double scale = 1.0;
    //to speed up rounding, scale=.1
    if (targetSum > 0){
        *(value) = -1*scale*(*(value));
    }
    else{
        *(value) = scale*(*(value));
    }
    //*(value)=5.0;
}



void copyGraph(double *A, double *Ain, mwSize n){
    for (int i = 0; i<n; i++){
        for (int j=0; j<n; j++){
            A[i*n+j]=Ain[i*n+j];
        }
    }

}



void copyIntoEmpty(int *list,int *result, int size){
    //result should be empty
    for (int i = 0; i < size; ++i){
        *(result+i) = *(list+i);
    }
}

int updateHashTable(value_edgeLocs_entry *valueLocPairs, double *value, int point){
    value_edgeLocs_entry *valueLocPair, *replaced;
    HASH_FIND(hh,valueLocPairs,value,sizeof(double),valueLocPair);
    int addFlag = 0;
    if (valueLocPair == NULL){
        valueLocPair = (value_edgeLocs_entry*)malloc(sizeof(value_edgeLocs_entry));
        valueLocPair->id = *(value);
        valueLocPair->edgeLocs = (int*)malloc(sizeof(int));
        valueLocPair->size = 1;
        valueLocPair->capacity = 1;
        valueLocPair->current = 0;
        *(valueLocPair->edgeLocs) = point;
        addFlag = 1;
        HASH_ADD(hh,valueLocPairs,id,sizeof(int),valueLocPair);                    
    }
    else{
        //add to the value loc pair
        if (valueLocPair->current+1 >= valueLocPair->size){
            int* toDelete = valueLocPair->edgeLocs;
            valueLocPair->edgeLocs = (int*)malloc(2*sizeof(int)*valueLocPair->size);
            //copy everything over 
            copyIntoEmpty(toDelete,valueLocPair->edgeLocs,valueLocPair->capacity);
            free(toDelete);
            valueLocPair->size = valueLocPair->size*2;
            valueLocPair->capacity = valueLocPair->capacity+1;
            valueLocPair->current = valueLocPair->current+1;
            *(valueLocPair->edgeLocs+valueLocPair->current)=point;
            HASH_REPLACE(hh,valueLocPairs,id,sizeof(int),valueLocPair,replaced);
        }
        else{
            valueLocPair->current = valueLocPair->current+1;
            valueLocPair->capacity = valueLocPair->capacity+1;
            *(valueLocPair->edgeLocs+valueLocPair->current)=point;
        }
    }
    return addFlag;
}

void round_forSpec(double *A, mwSize n, double *Ain, double targetSum, double *u, double *v, double *degree, double budget, double *amountRemovedInitialPr, double *threshPr, double *amountRemovedOutput, double *amountRemovedTotalPr, double *amountRemovedAfter, double *amountAddedPr){
    int i,j,k,l;
    int addFlag=0;
    int newFlag=0;
    int point = 0;
    double* min_key_add;
    double* min_key_sub;
    int* min_loc_add;
    int* min_loc_sub;
    value_edgeLocs_entry *valueLocPair, *replaced, *valueLocPairs=NULL;
    int m = (n*(n-1))/2;
    //storing where each of the edges is stored in data
    int* locations = (int*)malloc(m*sizeof(int));
    //two value for every edge, thes tart and end vertex
    int* data = (int*)malloc(m*2*sizeof(int));
    heap hAdd;
    heap_create(&hAdd,0,NULL);
    heap hSub;
    heap_create(&hSub,0,NULL);

    bool hashEmpty = 1;
    double thresh = *(threshPr);

    int numEdgesAdded = 0;
    int numEdgesRemoved = 0;


    copyGraph(A,Ain,n);
    double eigNorm = 0;

    for (i=0; i<n;++i){
        eigNorm = eigNorm + (u[i]*v[i]);
    }
    //ALL initialization
    for (i=0;i<n-1;++i){
        for (j=i+1;j<n;++j){
            //mexPrintf("edge %d %d \n",i,j);
            //for every pair, either they have the same sign or not
            double* value = (double*)malloc(sizeof(double));
            if (point > (m-1)){
                mexPrintf("point is too large \n");
                break;
            }
            //mexPrintf("point m %d %d \n", point,m);
            *(locations+point) = point;
            *(data+(2*point))=i;
            *(data+(2*point)+1)=j;
            //edge is on either side of the cut. 
            if ((u[i] < 0 && v[j]<0) || (u[i] > 0 && v[j]>0)){
                //if i need to sparisfy, add edges on either side of the cut. else, remove.
                if (targetSum < 0){
                    addFlag = 1;
                }
                compValueSingleEdge(value, u, v, degree, i, j, A, n, addFlag, eigNorm,targetSum); 
            }
            else{
                //if i need to strengthen cut, add edges crossing the cut. else, remove
                if (targetSum > 0){
                    addFlag = 1;
                }
                compValueSingleEdge(value, u, v, degree, i, j, A, n, addFlag, eigNorm,targetSum);
            }
            if (value == NULL){
                mexPrintf("value is null \n");
                break;
            }
            double absValue;
            if (*(value) <= 0){
                absValue = -1*(*(value));
            }
            else{
                absValue = *(value);
            }
            // code for saving space, if you only want to add edges to heap that have a significant 
            // edge change because those are the ones that will change the spectrum
            //int edgeChange = 0;
            //if (addFlag == 1 && ((1-A[i*n+j]) > .00001)){
            //    edgeChange = 1;
            //}
            //else if (addFlag ==0 && (A[i*n+j] > .00001)){
            //    edgeChange = 1;
            //}
            //if (absValue > 0.0000000000001 && edgeChange == 1){
            if (absValue > 0.0000000000001)
            // //lets us know if the value was new or not
                valueLocPair = NULL;
                if(hashEmpty == 1){
                    hashEmpty = 0;
                }
                HASH_FIND(hh,valueLocPairs,value,sizeof(double),valueLocPair);
                if (valueLocPair == NULL){
                    valueLocPair = (value_edgeLocs_entry*)malloc(sizeof(value_edgeLocs_entry));
                    if (value == NULL){
                        mexPrintf("value is null \n");
                        break;
                    }
                    valueLocPair->id = *(value);
                    //mexPrintf("value %f \n", *(value));
                    valueLocPair->edgeLocs = (int*)malloc(sizeof(int));
                    valueLocPair->size = 1;
                    valueLocPair->capacity = 1;
                    valueLocPair->current = 0;
                    if (valueLocPair->edgeLocs == NULL){
                        mexPrintf("edgelocs is null \n");
                    }
                    //else{
                        //mexPrintf("this is what edgelocs has in it %d \n",*(valueLocPair->edgeLocs));
                    //}
                    *(valueLocPair->edgeLocs) = point;
                    //mexPrintf("assigning point to edgelocs %d \n",*(valueLocPair->edgeLocs));
                    newFlag = 1;
                    HASH_ADD(hh,valueLocPairs,id,sizeof(double),valueLocPair); 
                    valueLocPair = NULL; 
                    HASH_FIND(hh,valueLocPairs,value,sizeof(double),valueLocPair);
                    if (valueLocPair == NULL){
                        mexPrintf("did not find what i just added \n");
                    }                  
                }
                else{
                //     //add to the value loc pair
                     if (valueLocPair->current+1 >= valueLocPair->size){
                        int* toDelete = valueLocPair->edgeLocs;
                        valueLocPair->edgeLocs = (int*)malloc(2*sizeof(int)*valueLocPair->size);
                        //copy everything over 
                        
                        copyIntoEmpty(toDelete,valueLocPair->edgeLocs,valueLocPair->capacity);
                        free(toDelete);
                        valueLocPair->size = valueLocPair->size*2;
                        valueLocPair->current = valueLocPair->current+1;
                        valueLocPair->capacity = valueLocPair->capacity+1;
                         *(valueLocPair->edgeLocs+valueLocPair->current)=point;
                         HASH_REPLACE(hh,valueLocPairs,id,sizeof(double),valueLocPair,replaced);
                     }
                     else{
                         valueLocPair->current = valueLocPair->current+1;
                         valueLocPair->capacity = valueLocPair->capacity+1;
                         *(valueLocPair->edgeLocs+valueLocPair->current)=point;
                     }
                 }
            }
            if (newFlag){
                if (addFlag){
                    //mexPrintf("new value in add! %f \n", *(value));
                    heap_insert(&hAdd, value, locations+point);
                }
                else{
                    //mexPrintf("new value in sub! %f \n", *(value));
                    heap_insert(&hSub, value, locations+point);
                }
            }
            newFlag = 0;
            addFlag = 0;
            point = point+1;
        }
    }
    //mexPrintf("exits for loop \n");
    // for (int k = 0; k< m; k++){
    //     mexPrintf("k %d \n", k);
    //     mexPrintf("i %d \n", (*(data+(2*k))));
    //     mexPrintf("j %d \n", (*(data+(2*k)+1)));
    // }


    //now i pop offf until value is at least target sum
    double eigDecrement = 0;
    double eigDecrement_add, eigDecrement_sub;
    double add1,add2,remove1,remove2,alpha;
    double amountAdded = 0;
    
    int count =0;
    int removing;
    double amountRemoved = *(amountRemovedTotalPr);

    if (amountRemoved>0){
        removing = 0;
    }
    else{
        removing = 1;
    }
    double amountRemovedAcross = 0;
    bool sub = heap_delmin(&hSub, (void**)&min_key_sub, (void**)&min_loc_sub);
    count = 0;

    bool add = heap_delmin(&hAdd, (void**)&min_key_add, (void**)&min_loc_add);
    count = 0;
    
    //mexPrintf("target %f\n", targetSum);
    int itercount = 0;
    if (targetSum > 0){
        targetSum = -1*targetSum;
    }
    //delta - s < -thresh -> delta 
    while ((targetSum - eigDecrement < -thresh) && (sub && add) && (amountRemovedAcross < budget)) {
        if (itercount % 100 == 0){
            //mexPrintf("itercount %d \n", itercount);
            //mexPrintf("target - eigDecrment Sum %f \n", targetSum-eigDecrement);
            //smexPrintf("budget - amountRemovedAcross %f \n", budget-amountRemovedAcross);

        }
        itercount = itercount+1;
    //while (count < 6000){
        
        //remove things until i 
        if (removing == 1){
                //mexPrintf("removing \n");
                if (min_key_sub==NULL){
                    mexPrintf("min key sub \n");
                    break;
                }

                
                eigDecrement_sub = (*(min_key_sub));
                
                //mexPrintf("amount removed, budget %f %f\n", amountRemoved, budget);
                if (eigDecrement_sub >= -0.000000001){
                    mexPrintf("hit zero valued keys\n");
                    break;
                }
                //making sure that current score + new score is not more negative than target sum, 
                //could be all of them break it
                if (eigDecrement_sub + eigDecrement >= targetSum){
                    //mexPrintf("removing \n");
                    //mexPrintf("eigDecrement sub, target sum, eigDecrement %f %f %f\n", eigDecrement_sub, targetSum, eigDecrement);
                    valueLocPair = NULL;
                    HASH_FIND(hh,valueLocPairs,min_key_sub,sizeof(double),valueLocPair);
                    if (valueLocPair == NULL){
                        mexPrintf("value is null \n");
                        break;
                    }
                    int* locPr = valueLocPair->edgeLocs+valueLocPair->current;
                    if (locPr == NULL){
                        mexPrintf("loc pr is null \n");
                        break;
                    }
                    int loc = *(locPr);
                    
                    if (loc > m-1){
                        mexPrintf("loc large \n");
                        break;
                    }
                    eigDecrement = eigDecrement + eigDecrement_sub;
                    k = *(data+(2*loc));
                    l = *(data+(2*loc)+1);
                    if (k*n+l > (n*n)-1){
                        mexPrintf("vertex not real \n");
                        break;
                    }
                    if (l*n+k > (n*n)-1){
                        mexPrintf("vertex not real \n");
                        break;
                    }
                    sub = (heap_delmin(&hSub, (void**)&min_key_sub, (void**)&min_loc_sub)); 
                    
                    //make sure this was not the last edge. so im checking, if it was the last one
                    //i do NOT remove it. 
                    
                    if (sub){
                        double prev = A[k*n+l];
                        amountRemoved = amountRemoved + prev;
                        amountRemovedAcross = amountRemovedAcross + prev;
                        //mexPrintf("amount removed %f\n", A[k*n+l]);
                        //mexPrintf("coordinates of edge removed %d %d \n",k,l);
                        A[k*n+l]=0;
                        A[l*n+k]=0;
                        numEdgesRemoved = numEdgesRemoved+1;
                        //mexPrintf("edge degrees eigvecEntries graphBefore %d %d %f %f %f %f %f",k,l,degree[k],degree[l],u[k],v[l],prev);
                        if (amountRemoved > 0){
                            removing = 0;
                            if (valueLocPair->current > 0){
                                //get next edge with this value
                                valueLocPair->current = valueLocPair->current-1;
                                //put this back on the heap
                                heap_insert(&hSub, &eigDecrement_sub, 0);
                            } 
                        }
                    }
                    else{
                        break;
                    }
                         
                }
                else{
                    //means i was too high. still need to replace the edge. also break if sub isn't valid
                    sub = (heap_delmin(&hSub, (void**)&min_key_sub, (void**)&min_loc_sub));
                    if(!sub){
                        break;
                    }
                }
                //sub = (heap_delmin(&hSub, (void**)&min_key_sub, (void**)&min_loc_sub));
        }
        else {
            //mexPrintf("adding \n");
            if (min_key_add == NULL){
                mexPrintf("add null \n");
                break;
            }
            eigDecrement_add = (*(min_key_add));
            if (eigDecrement_add + eigDecrement >= targetSum){
                //mexPrintf("eig decrement add, eigDecrement, targetSum %f %f %f \n", eigDecrement_add, eigDecrement, targetSum);
                if (eigDecrement_add >= -0.000000001){
                    mexPrintf("zero\n");
                    break;
                }
                valueLocPair = NULL;
                HASH_FIND(hh,valueLocPairs,min_key_add,sizeof(double),valueLocPair);
                if (valueLocPair == NULL){
                    mexPrintf("value is null \n");
                    break;
                }
                int* locPr = valueLocPair->edgeLocs+valueLocPair->current;
                if (locPr == NULL){
                    mexPrintf("loc pr is null \n");
                    break;
                }
                int loc = *(locPr);

                if (loc > m-1){
                    mexPrintf("loc null \n");
                    break;
                }
                if (eigDecrement_add >= 0){
                        mexPrintf("zero\n");
                        break;
                    }
                eigDecrement = eigDecrement + eigDecrement_add;
                i = *(data+(2*loc));
                j = *(data+(2*loc)+1);
                if (i*n+j > (n*n)-1){
                    mexPrintf("vertex not real \n");
                    break;
                }
                if (j*n+i > (n*n)-1){
                    mexPrintf("vertex not real \n");
                    break;
                }
                amountRemoved = amountRemoved - (1-A[i*n+j]);
                //mexPrintf("adding %d %d %f %f \n",i,j,u[i],u[j]);
                //mexPrintf("amount removed %f\n", amountRemoved);
                //mexPrintf("amount added %f\n", 1-A[i*n+j]);
                amountAdded = amountAdded + (1-A[i*n+j]);
                A[i*n+j]=1;
                A[j*n+i]=1;
                numEdgesAdded = numEdgesAdded + 1;
                //mexPrintf("edge degrees eigvecEntries amountAdded %d %d %f %f %f %f %f",i,j,degree[i],degree[j],u[i],v[j],(1-A[i*n+j]));
                if (amountRemoved <=0){
                    removing = 1;
                    if (valueLocPair->current > 0){
                            //get next edge with this value
                            valueLocPair->current = valueLocPair->current-1;
                            //put this back on the heap
                            heap_insert(&hAdd, &eigDecrement_add, 0);
                    } 
                }
            }
            else{
                //mexPrintf("could not add \n");
            }
            add = heap_delmin(&hAdd, (void**)&min_key_add, (void**)&min_loc_add);
        }
        count++;
    }
   
    *(amountRemovedOutput) = amountRemovedAcross;
    *(amountRemovedAfter) = amountRemoved;
    *(amountAddedPr) = amountAdded;
    heap_destroy(&hAdd);
    heap_destroy(&hSub);
    struct value_edgeLocs_entry *current_user, *tmp;

    HASH_ITER(hh, valueLocPairs, current_user, tmp) {
         HASH_DEL(valueLocPairs,current_user);  /* delete; users advances to next */
         free(current_user);            /* optional- if you want to free  */
    }

}


    

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *out;
    mxArray *amountRemovedOutput;
    mxArray *amountRemovedTotalOutput;
    mxArray *amountAddedOutput;
    mxArray *removeFlag;
    double *Ain;
    double *targetSumArray;
    double *budgetArray;
    double targetSum;
    double budget;
    double thresh;

    double *u;
    double *v;
    double *degree;

    double *A;
    double *amountRemovedPr;
    double *amountRemovedInitialPr;
    double *amountAddedPr;
    double *threshPr;

    double *amountRemovedTotalPr;
    double *amountRemovedAfter;
    double *amountRemovedAcrossPr;

    
    
    Ain= mxGetPr(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    targetSumArray = mxGetPr(prhs[1]);
    u = mxGetPr(prhs[2]);
    v = mxGetPr(prhs[3]);
    degree = mxGetPr(prhs[4]);
    budgetArray = mxGetPr(prhs[5]);
    amountRemovedAcrossPr = mxGetPr(prhs[6]);
    threshPr = mxGetPr(prhs[7]);
    amountRemovedTotalPr = mxGetPr(prhs[8]);
    targetSum = targetSumArray[0];
    budget = budgetArray[0];

    out = plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
    amountRemovedOutput = plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    amountRemovedTotalOutput = plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    amountAddedOutput = plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);

    A = mxGetPr(out);
    amountRemovedPr = mxGetPr(amountRemovedOutput);
    amountRemovedAfter = mxGetPr(amountRemovedTotalOutput);
    amountAddedPr = mxGetPr(amountAddedOutput);




    //mexPrintf("%d\n", n);

    //double alpha = compAlpha(0.0,0.0,1.0,1.0);
    //mexPrintf("alpha testing \n");

    //genEqualityConstraints(A,n,Ain,targetSum,u,v,degree);
    round_forSpec(A,n,Ain,targetSum,u,v,degree,budget,amountRemovedAcrossPr,threshPr,amountRemovedPr,amountRemovedTotalPr,amountRemovedAfter,amountAddedPr);


}
