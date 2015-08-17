#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;

// line for multiple processing: #pragma parallel =>  OpenMP Cheat Sheet
  
//TODO: tidy up functions - within 80char/line
  
  
//####################################################################
//####################################################################
//####################################################################

//####################################################################
// [[Rcpp::export]]
NumericVector inertiaCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> sender,
  std::vector<std::string> target,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string inertiatype ) {
  
  NumericVector result(sender.size());
  
  // for-loop i: for each event, do:
    for ( int i = 0; i < sender.size(); i++){
      
      // reset all the variables
      double inertia = 0;
      double totalinertia = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // sender-target inertia only
        if ( inertiatype == "s-t-only" ){
          if ( target[w] == target[i] && sender[w] == sender[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch 
        if ( inertiatype == "s-t-typematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix
        if ( inertiatype == "s-t-typemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter
        if ( inertiatype == "s-t-typefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributematch
        if ( inertiatype == "s-t-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributemix
        if ( inertiatype == "s-t-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributefilter
        if ( inertiatype == "s-t-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributematch
        if ( inertiatype == "s-t-typematch-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributemix
        if ( inertiatype == "s-t-typematch-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributefilter
        if ( inertiatype == "s-t-typematch-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributematch
        if ( inertiatype == "s-t-typemix-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributemix
        if ( inertiatype == "s-t-typemix-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributefilter
        if ( inertiatype == "s-t-typemix-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributematch
        if ( inertiatype == "s-t-typefilter-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributemix
        if ( inertiatype == "s-t-typefilter-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributefilter
        if ( inertiatype == "s-t-typefilter-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totalinertia = totalinertia + 0;
        } else {
          inertia = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totalinertia = totalinertia + inertia;
        }
        
        resulttemp = totalinertia;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable		
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}


//####################################################################
// [[Rcpp::export]]
NumericVector degreeCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> degreevar,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string degreetype ) {
  
  NumericVector result(degreevar.size());
  
  // for-loop i: for each event, do:
    for ( int i = 0; i < degreevar.size(); i++){
      
      // reset all the variables
      double degree = 0;
      double totaldegree = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // degreevar only
        if ( degreetype == "d-only" ){
          if ( degreevar[w] == degreevar[i]  ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch 
        if ( degreetype == "d-typematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix
        if ( degreetype == "d-typemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter
        if ( degreetype == "d-typefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributematch
        if ( degreetype == "d-attributematch" ){
          if ( degreevar[w] == degreevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributemix
        if ( degreetype == "d-attributemix" ){
          if ( degreevar[w] == degreevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributefilter
        if ( degreetype == "d-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributematch
        if ( degreetype == "d-typematch-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributemix
        if ( degreetype == "d-typematch-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributefilter
        if ( degreetype == "d-typematch-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributematch
        if ( degreetype == "d-typemix-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributemix
        if ( degreetype == "d-typemix-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributefilter
        if ( degreetype == "d-typemix-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributematch
        if ( degreetype == "d-typefilter-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributemix
        if ( degreetype == "d-typefilter-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributefilter
        if ( degreetype == "d-typefilter-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totaldegree = totaldegree + 0;
        } else {
          degree = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totaldegree = totaldegree + degree;
        }
        
        resulttemp = totaldegree;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable		
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}



//####################################################################
// [[Rcpp::export]]
NumericVector fourCycleCpp(
  std::vector<std::string> sender,
  std::vector<std::string> target,
  std::vector<std::string> typevar,
  NumericVector time,
  NumericVector weightvar,
  double xlog,
  std::vector<std::string> attrvarNow,
  std::string attrNow,
  std::vector<std::string> attrvarAaj,
  std::string attrAaj, 
  std::vector<std::string> attrvarBib,
  std::string attrBib,
  std::vector<std::string> attrvarCij,
  std::string attrCij, 
  std::string fourCycleType) {
  
  NumericVector result(sender.size());
  std::vector<std::string> w;
  std::vector<std::string> x;
  std::vector<std::string> y;
  std::vector<std::string> wy;
  double weightA;
  double weightB;
  double weightC;
  double tempTotalWeightC;
  double tempTotalWeightA;
  double tempTotalWeightAPositive;
  double tempTotalWeightANegative;
  double tempTotalWeightB;
  double tempTotalWeightCPositive;
  double tempTotalWeightCNegative;
  double tempTotalWeightABC;
  double tempTotalWeightABCPositive;
  double tempTotalWeightABCNegative;
  double totalWeightABCPositive;
  double totalWeightABCNegative;
  double totalWeightABC;
  
  
  //for each event (=i-loop open)
  for ( int i = 0; i < sender.size(); i++){
    
    //Filter: only current events with given attribute are selected.
    if ( attrvarNow[i] == attrNow ) {
      
      w.clear();
      x.clear();
      y.clear();
      totalWeightABCPositive = 0;
      totalWeightABCNegative = 0;
      tempTotalWeightABC = 0;
      tempTotalWeightABCNegative = 0;
      tempTotalWeightABCPositive = 0;
      totalWeightABC = 0;
      
      // generate a list of issues that $a$ has used in past (j-loop-open) (=w-vector)
      for ( int j = 0; j < i; j++ ) {
        if (sender[j] == sender[i] && target[j] != target[i] && attrvarAaj[j] == attrAaj){
          w.push_back(target[j]);
        }
        // list of actors who also used p (the same type/way a used it)
        // here: attrvarM1[j] == attr2 == attr1 für nodematch
        // here: attrvarM1[j] == attr2 != attr1 für nodemix
        if (fourCycleType == "standard"){
          if (sender[j] != sender[i] && target[j] == target[i] && attrvarBib[j] == attrBib ){
            x.push_back(sender[j]);
          }
        }//closes if-fourCycleType == standard
        if (fourCycleType == "positive"){
          if (sender[j] != sender[i] && target[j] == target[i] && typevar[j] == typevar[i] && attrvarBib[j] == attrBib ){
            x.push_back(sender[j]);
          }
        }//closes if-fourCycleType == positive (=supporting)
        if (fourCycleType == "negative"){
          if (sender[j] != sender[i] && target[j] == target[i] && typevar[j] != typevar[i] && attrvarBib[j] == attrBib ){
            x.push_back(sender[j]);
          }
        }//closes if-fourCycleType == negative (=opposing)
      } // j-loop close
      
      // clean up w (only unique values)
      std::sort( w.begin(), w.end() );
      w.erase( unique( w.begin(), w.end() ), w.end() );
      // clean up x (only unique values)
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each person in the list x (m-loop open) (=y-vector; wy-vector)
      for (int m = 0; m < x.size(); m++ ) {
        
        tempTotalWeightB = 0;
        
        // What did actor say in past? (n-loop open)
        y.clear();
        for ( int n = 0; n < i; n++ ) {              
          // for each person: find y-vector (list of targets $i$ has used)
          if (sender[n] != sender[i] && sender[n] == x[m] && target[n] != target[i] && attrvarCij[n] == attrCij ){
            y.push_back(target[n]);
          }
        } // n-loop close
        
        // clean up y (only unique values)
        std::sort( y.begin(), y.end() );
        y.erase( unique( y.begin(), y.end() ), y.end() );
        
        // interlock between x and y = wy vector
        wy.clear();
        std::sort( w.begin(), w.end() );
        std::sort( y.begin(), y.end() );
        std::set_intersection (w.begin(), w.end(), y.begin(), y.end(), std::back_inserter(wy) );
        // erase douplicates
        sort( wy.begin(), wy.end() );
        wy.erase( unique( wy.begin(), wy.end() ), wy.end() );
        
        // if these two actors (x-vector from a and y-vector from i) interact via a shared concept: calculate weightB (weightB = w_t(i, p))
        // positive/negative four cycle: only choose those events with different type
        if ( fourCycleType == "standard" ){
          if (wy.size() != 0) {
            for ( int o = 0; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == target[i] && attrvarBib[o] == attrBib ) {
                weightB = std::abs(weightvar[o]) * exp( - ( time[i] - time[o] ) * xlog)  * xlog;
                if ( time[i] == time[o] ) {
                  weightB = 0;
                }
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == standard
        if ( fourCycleType == "positive" ){
          if (wy.size() != 0) {
            for ( int o = 0; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == target[i] && typevar[o] == typevar[i] && attrvarBib[o] == attrBib ) {
                weightB = std::abs(weightvar[o]) * exp( - ( time[i] - time[o] ) * xlog)  * xlog;
                if ( time[i] == time[o] ) {
                  weightB = 0;
                }
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == positive
        if ( fourCycleType == "negative" ){
          if (wy.size() != 0) {
            for ( int o = 0; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == target[i] && typevar[o] != typevar[i] && attrvarBib[o] == attrBib ) {
                weightB = std::abs(weightvar[o]) * exp( - ( time[i] - time[o] ) * xlog)  * xlog;
                if ( time[i] == time[o] ) {
                  weightB = 0;
                }
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == negative
        
        if ( fourCycleType == "standard"){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightA = 0;
            tempTotalWeightC = 0;
            
            for ( int q = 0; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && attrvarCij[q] == attrCij ) {
                weightC = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightC = 0;
                }
                tempTotalWeightC = tempTotalWeightC + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == sender[i] && target[q] == wy[p] && attrvarAaj[q] == attrAaj ) {
                weightA = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightA = 0;
                }
                tempTotalWeightA = tempTotalWeightA + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABC = tempTotalWeightA * tempTotalWeightB * tempTotalWeightC;
            
            // for each person: for each entry = sum up the multiplications
            totalWeightABC = tempTotalWeightABC + totalWeightABC;
          }//closes p-loop
        }//closes if fourCycleType == standard
        if ( fourCycleType == "positive" ){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightAPositive = 0;
            tempTotalWeightCPositive = 0;
            tempTotalWeightANegative = 0;
            tempTotalWeightCNegative = 0;
            
            for ( int q = 0; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // if both weightA and weightC events are of same type, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] == typevar[i] && attrvarCij[q] == attrCij ) {
                weightC = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightC = 0;
                }
                tempTotalWeightCPositive = tempTotalWeightCPositive + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == sender[i] && target[q] == wy[p] && typevar[q] == typevar[i] && attrvarAaj[q] == attrAaj ) {
                weightA = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightA = 0;
                }
                tempTotalWeightAPositive = tempTotalWeightAPositive + weightA;
              }
              // if both weightA and weightC events are of same type and negative, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] != typevar[i] && attrvarCij[q] == attrCij ) {
                weightC = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightC = 0;
                }
                tempTotalWeightCNegative = tempTotalWeightCNegative + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == sender[i] && target[q] == wy[p] && typevar[q] != typevar[i] && attrvarAaj[q] == attrAaj ) {
                weightA = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog) * xlog;
                if ( time[i] == time[q] ) {
                  weightA = 0;
                }
                tempTotalWeightANegative = tempTotalWeightANegative + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABCPositive = tempTotalWeightAPositive * tempTotalWeightB * tempTotalWeightCPositive;
            tempTotalWeightABCNegative = tempTotalWeightANegative * tempTotalWeightB * tempTotalWeightCNegative;
            
            // for each person: for each entry = sum up the multiplications
            totalWeightABCPositive = tempTotalWeightABCPositive + totalWeightABCPositive;
            totalWeightABCNegative = tempTotalWeightABCNegative + totalWeightABCNegative;
          }//closes p-loop
        }//closes if fourCycleType == positive
        if ( fourCycleType == "negative" ){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightAPositive = 0;
            tempTotalWeightCPositive = 0;
            tempTotalWeightANegative = 0;
            tempTotalWeightCNegative = 0;
            
            for ( int q = 0; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // if both weightA and weightC events are of opposite type, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] == typevar[i] && attrvarCij[q] == attrCij ) {
                weightC = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightC = 0;
                }
                tempTotalWeightCPositive = tempTotalWeightCPositive + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == sender[i] && target[q] == wy[p] && typevar[q] != typevar[i] && attrvarAaj[q] == attrAaj ) {
                weightA = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightA = 0;
                }
                tempTotalWeightAPositive = tempTotalWeightAPositive + weightA;
              }
              // if both weightA and weightC events are negative, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] != typevar[i] && attrvarCij[q] == attrCij ) {
                weightC = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog)  * xlog;
                if ( time[i] == time[q] ) {
                  weightC = 0;
                }
                tempTotalWeightCNegative = tempTotalWeightCNegative + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == sender[i] && target[q] == wy[p] && typevar[q] == typevar[i] && attrvarAaj[q] == attrAaj ) {
                weightA = std::abs(weightvar[q]) * exp( - ( time[i] - time[q] ) * xlog) * xlog;
                if ( time[i] == time[q] ) {
                  weightA = 0;
                }
                tempTotalWeightANegative = tempTotalWeightANegative + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABCPositive = tempTotalWeightAPositive * tempTotalWeightB * tempTotalWeightCPositive;
            tempTotalWeightABCNegative = tempTotalWeightANegative * tempTotalWeightB * tempTotalWeightCNegative;
            
            // for each person: for each entry = sum up the multiplications
            totalWeightABCPositive = tempTotalWeightABCPositive + totalWeightABCPositive;
            totalWeightABCNegative = tempTotalWeightABCNegative + totalWeightABCNegative;
          }//closes p-loop
        }//closes if fourCycleType == negative
      } // m-loop close
      
      if ( fourCycleType == "standard"){
        totalWeightABC = std::pow(totalWeightABC, 1/3.); //whyever, there has to be a . behind the 1/3
      }else{
        totalWeightABC = std::pow(totalWeightABCPositive, 1/3.) + std::pow(totalWeightABCNegative, 1/3.); //wieso auch immer - aber da muss ein Punkt hinter die 3
      }
      result[i] = totalWeightABC; 
      
    }else{ //closes "if ( attrvarNow[i] == attrNow ) {}"
           result[i] = 0;
    }
  } // i-loop close
  return Rcpp::wrap(result);
}


//####################################################################
// [[Rcpp::export]]
NumericVector similarityTotalAverageCpp(
  std::vector<std::string> sender,
  std::vector<std::string> target,
  NumericVector time,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string totalAverageSim,
  std::string matchNomatchSim,
  std::string senderTargetSim) {
  
  NumericVector result(sender.size());
  std::vector<std::string> v;
  std::vector<std::string> w;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  std::vector<double> a_positive;
  std::vector<double> a_negative;
  std::vector<double> i_negative;
  std::vector<double> i_positive;	
  std::vector<std::string> xwneg;
  double totalNumber;
  double numberNoMatch;
  double numberMatch;
  
  //for each event (i-loop open)
  for ( int i = 0; i < sender.size(); i++){
    
    v.clear();
    w.clear();
    xw.clear();
    i_negative.clear();
    a_negative.clear();
    i_positive.clear();
    a_positive.clear();
    xwneg.clear();
    totalNumber = 0;
    
    // loop back and find actors who said p and concepts that a said
    for (int j = 0; j < i; j++) {
      if ( senderTargetSim == "sender"){
        if ( matchNomatchSim == "match" ){
          if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute && eventTypeVar[j] == eventTypeVar[i]){
            v.push_back(sender[j]);   // v = who else has said p just like $a$ used p (typematch)?
          }
        }else{ //no match-option
               if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
                 v.push_back(sender[j]);   // v = who else has said p?
               }
        }
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(target[j]);   // w = what else has actor a said?
        }
        
      }else{
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          v.push_back(target[j]);   // v = what else has actor a said?
        }
        if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(sender[j]);   // w = who else has said p?
        }
      }
    }
    // clean up v and w
    std::sort( w.begin(), w.end() );
    w.erase( unique( w.begin(), w.end() ), w.end() );
    std::sort( v.begin(), v.end() );
    v.erase( unique( v.begin(), v.end() ), v.end() );
    
    // for each entry in v
    for (int k = 0; k < v.size(); k++){
      
      x.clear();
      xw.clear();
      
      if ( senderTargetSim == "sender" ){
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (sender[l] == v[k] && time[l] != time[i] && target[l] != target[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(target[l]);
          }
        }//closes l-loop
      }else{
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (target[l] == v[k] && time[l] != time[i] && sender[l] != sender[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(sender[l]);
          }
        }//closes l-loop
      }
      
      // clean up x-vector
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
      xw.clear();
      std::sort( w.begin(), w.end() );
      std::sort( x.begin(), x.end() );
      std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
      // erase douplicates
      sort( xw.begin(), xw.end() );
      xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
      
      if ( matchNomatchSim == "match"){
        
        // if there acctually is an intersection in xw
        if (xw.size() != 0 ) {
          
          // for each entry in xw:
            for (int m = 0; m < xw.size(); m++) {
              
              i_negative.clear();
              i_positive.clear();
              a_negative.clear();
              a_positive.clear();
              
              // loop back over all events until i
              for (int n = 0; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }//closes senderTargetSim == target
              }//closes n-loop		
            }//closes m-loop
          
        }//closes if xw.size != 0
        
        numberNoMatch = 0;
        numberMatch = 0;
        
        //how large are the respextive vectors with the matches in them?
        if (a_positive.size() >= i_positive.size() && i_positive.size() != 0 ) {
          numberMatch = i_positive.size();
        }
        if (a_negative.size() >= i_negative.size() && i_negative.size() != 0) {
          numberNoMatch = i_negative.size();
        }
        if (i_positive.size() > a_positive.size() && a_positive.size() != 0) {
          numberMatch = a_positive.size();
        }
        if (i_negative.size() > a_negative.size() && a_negative.size() != 0) {
          numberNoMatch = a_negative.size();
        }
        //// how many actors used concept v[k] in same manner as a? // how many concepts are used by both in same manner?
        totalNumber = totalNumber + numberNoMatch + numberMatch;
      }else{ // if matchNomatchSim = "nomatch"
             // how many actors used concept v[k]? // how many concepts are used by both?
             totalNumber = totalNumber + xw.size();
      }
      
    }//closes k-loop
    
    if (v.size() != 0){
      if (totalAverageSim == "total") {
        result[i] = totalNumber;
      }
      if (totalAverageSim == "average") {
        result[i] = totalNumber/v.size(); //TODO: correct to divide by v.size?
      }
    }else{
      result[i] = 0;
    }
  }//i-loop close
  return Rcpp::wrap(result);
}

//####################################################################
// [[Rcpp::export]]
NumericVector similaritySimpleCpp(
  std::vector<std::string> sender,
  std::vector<std::string> target,
  NumericVector time,
  double xlog,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string matchNomatchSim,
  std::string senderTargetSim) {
  
  NumericVector result(sender.size());
  std::vector<std::string> v;
  std::vector<std::string> w;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  std::vector<double> a_positive;
  std::vector<double> a_negative;
  std::vector<double> i_negative;
  std::vector<double> i_positive;	
  std::vector<std::string> xwneg;
  double totalNumber;
  double numberNoMatch;
  double numberMatch;
  double timePLast;
  int counter;
  double totalSim;
  double weightSim;
  
  //for each event (i-loop open)
  for ( int i = 0; i < sender.size(); i++){
    
    v.clear();
    w.clear();
    xw.clear();
    i_negative.clear();
    a_negative.clear();
    i_positive.clear();
    a_positive.clear();
    xwneg.clear();
    totalNumber = 0;
    
    // loop back and find actors who said p and concepts that a said
    for (int j = 0; j < i; j++) {
      if ( senderTargetSim == "sender"){
        if ( matchNomatchSim == "match" ){
          if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute && eventTypeVar[j] == eventTypeVar[i]){
            v.push_back(sender[j]);   // v = who else has said p just like $a$ used p (typematch)?
          }
        }else{ //no match-option
               if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
                 v.push_back(sender[j]);   // v = who else has said p?
               }
        }
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(target[j]);   // w = what else has actor a said?
        }
        
      }else{
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          v.push_back(target[j]);   // v = what else has actor a said?
        }
        if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(sender[j]);   // w = who else has said p?
        }
      }
    }
    // clean up v and w
    std::sort( w.begin(), w.end() );
    w.erase( unique( w.begin(), w.end() ), w.end() );
    std::sort( v.begin(), v.end() );
    v.erase( unique( v.begin(), v.end() ), v.end() );
    
    // set variables
    totalSim = 0;
    counter = 0;
    
    // for each entry in v
    for (int k = 0; k < v.size(); k++){
      
      x.clear();
      xw.clear();
      
      if ( senderTargetSim == "sender" ){
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (sender[l] == v[k] && time[l] != time[i] && target[l] != target[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(target[l]);
          }
        }//closes l-loop
      }else{
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (target[l] == v[k] && time[l] != time[i] && sender[l] != sender[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(sender[l]);
          }
        }//closes l-loop
      }
      
      // clean up x-vector
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
      xw.clear();
      std::sort( w.begin(), w.end() );
      std::sort( x.begin(), x.end() );
      std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
      // erase douplicates
      sort( xw.begin(), xw.end() );
      xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
      
      if ( matchNomatchSim == "match"){
        
        // if there acctually is an intersection in xw
        if (xw.size() != 0 ) {
          
          // for each entry in xw:
            for (int m = 0; m < xw.size(); m++) {
              
              i_negative.clear();
              i_positive.clear();
              a_negative.clear();
              a_positive.clear();
              
              // loop back over all events until i
              for (int n = 0; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }//closes senderTargetSim == target
              }//closes n-loop		
            }//closes m-loop
          
        }//closes if xw.size != 0
        
        numberNoMatch = 0;
        numberMatch = 0;
        
        //how large are the respextive vectors with the matches in them?
        if (a_positive.size() >= i_positive.size() && i_positive.size() != 0 ) {
          numberMatch = i_positive.size();
        }
        if (a_negative.size() >= i_negative.size() && i_negative.size() != 0) {
          numberNoMatch = i_negative.size();
        }
        if (i_positive.size() > a_positive.size() && a_positive.size() != 0) {
          numberMatch = a_positive.size();
        }
        if (i_negative.size() > a_negative.size() && a_negative.size() != 0) {
          numberNoMatch = a_negative.size();
        }
        //// how many actors used concept v[k] in same manner as a? // how many concepts are used by both in same manner?
        totalNumber = totalNumber + numberNoMatch + numberMatch;
      }else{ // if matchNomatchSim = "nomatch"
             // how many actors used concept v[k]? // how many concepts are used by both?
             totalNumber = totalNumber + xw.size();
      }//closes matchNomatchSim == "nomatch"
      
      // how many actors used concept v[k]? // how many concepts are used by both?
      if ( xw.size() != 0 ){
        // find time for the time-discount in the simple-similarity equation
        if ( senderTargetSim == "sender"){
          // find time, when actor k used concept $p$ last
          for (int q = i-1; q >= 0; q--){
            if ( matchNomatchSim == "nomatch"){
              if (sender[q] == v[k] && target[q] == target[i] && time[q] != time[i] &&
                    eventAttributeVar[q] == eventAttribute){
                timePLast = time[q];
                break;
              }
            }else if (matchNomatchSim == "match"){
              if (sender[q] == v[k] && target[q] == target[i] && time[q] != time[i] &&
                    eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == eventTypeVar[i]){
                timePLast = time[q];
                break;
              }
            }
          }//closes q-loop
        }else if ( senderTargetSim == "target"){
          // find time, when target k was last used by $a$
            for (int q = i-1; q >= 0; q--){
              if (matchNomatchSim == "nomatch"){
                if (target[q] == v[k] && sender[q] == sender[i] && time[q] != time[i] &&
                      eventAttributeVar[q] == eventAttribute){
                  timePLast = time[q];
                  break;
                }
              } else if (matchNomatchSim == "match"){
                if (target[q] == v[k] && sender[q] == sender[i] && time[q] != time[i] &&
                      eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == eventTypeVar[i]){
                  timePLast = time[q];
                  break;
                }
              }
            }//closes q-loop
        }//closes senderTargetSim == "target"
        
        // calculate weight of each actor/target
        weightSim = totalNumber * exp(-(time[i]-timePLast)*xlog) * xlog;
        totalSim = totalSim + weightSim;
        if (weightSim != 0) {
          counter++;   
        }
      }//closes if xw.size != 0
    }//closes k-loop
    
    if (counter == 0) {
      result[i] = 0;
    }else{
      result[i] = (totalSim/counter);
    }
  }//i-loop close
  return Rcpp::wrap(result);
}

//####################################################################
// [[Rcpp::export]]
NumericVector similarityComplexCpp(
  std::vector<std::string> sender,
  std::vector<std::string> target,
  NumericVector time,
  double xlog,
  double halflifeTimeDifference,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string matchNomatchSim,
  std::string senderTargetSim) {
  
  NumericVector result(sender.size());
  std::vector<std::string> v;
  std::vector<std::string> w;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  std::vector<std::string> xwneg;
  std::vector<double> a_positive;
  std::vector<double> a_negative;
  std::vector<double> i_negative;
  std::vector<double> i_positive;	
  std::vector<double> i_sendertarget;
  std::vector<double> a_sendertarget;
  double totalNumber;
  double timePLast;
  int counter;
  double totalSim;
  double weightSim;
  double sumCouplePositive;
  double sumCoupleNegative;
  double couple;
  double sumCoupleNoMatch;
  double sumConcept;
  
  
  //for each event (i-loop open)
  for ( int i = 0; i < sender.size(); i++){
    
    v.clear();
    w.clear();
    xw.clear();
    i_negative.clear();
    a_negative.clear();
    i_positive.clear();
    a_positive.clear();
    a_sendertarget.clear();
    i_sendertarget.clear();
    xwneg.clear();
    totalNumber = 0;
    counter = 0;
    sumConcept = 0; //TODO: here or inside k-loop?
    
    
    // loop back and find actors who said p and concepts that a said
    for (int j = 0; j < i; j++) {
      if ( senderTargetSim == "sender"){
        if ( matchNomatchSim == "match" ){
          if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute && eventTypeVar[j] == eventTypeVar[i]){
            v.push_back(sender[j]);   // v = who else has said p just like $a$ used p (typematch)?
          }
        }else{ //no match-option
               if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
                 v.push_back(sender[j]);   // v = who else has said p?
               }
        }
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(target[j]);   // w = what else has actor a said?
        }
        
      }else{
        if (sender[j] == sender[i] && target[j] != target[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          v.push_back(target[j]);   // v = what else has actor a said?
        }
        if (target[j] == target[i] && sender[j] != sender[i] && time[j] != time[i] && eventAttributeVar[j] == eventAttribute){
          w.push_back(sender[j]);   // w = who else has said p?
        }
      }
    }
    // clean up v and w
    std::sort( w.begin(), w.end() );
    w.erase( unique( w.begin(), w.end() ), w.end() );
    std::sort( v.begin(), v.end() );
    v.erase( unique( v.begin(), v.end() ), v.end() );
    
    // set variables
    totalSim = 0;
    counter = 0;
    
    // for each entry in v
    for (int k = 0; k < v.size(); k++){
      
      x.clear();
      xw.clear();
      
      if ( senderTargetSim == "sender" ){
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (sender[l] == v[k] && time[l] != time[i] && target[l] != target[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(target[l]);
          }
        }//closes l-loop
      }else{
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (target[l] == v[k] && time[l] != time[i] && sender[l] != sender[i] && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(sender[l]);
          }
        }//closes l-loop
      }
      
      // clean up x-vector
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
      xw.clear();
      std::sort( w.begin(), w.end() );
      std::sort( x.begin(), x.end() );
      std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
      // erase douplicates
      sort( xw.begin(), xw.end() );
      xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
      
      if ( matchNomatchSim == "match"){
        
        // if there acctually is an intersection in xw
        if (xw.size() != 0 ) {
          
          // for each entry in xw:
            for (int m = 0; m < xw.size(); m++) {
              
              i_negative.clear();
              i_positive.clear();
              a_negative.clear();
              a_positive.clear();
              i_sendertarget.clear();
              a_sendertarget.clear();
              
              // loop back over all events until i
              for (int n = 0; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (sender[n] == sender[i] && target[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] == eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (target[n] == target[i] && sender[n] == xw[m] && eventTypeVar[n] != eventTypeVar[i] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }//closes senderTargetSim == target
              }//closes n-loop		
            }//closes m-loop
          
        }//closes if xw.size != 0
        
        sumCouplePositive = 0;
        sumCoupleNegative = 0;
        couple = 0;
        
        //how large are the respextive vectors with the matches in them?
        if (a_positive.size() >= i_positive.size() && i_positive.size() != 0 ) {
          for (int p = 0; p < i_positive.size(); p++){
            couple = 1 * exp(-(std::abs(i_positive[p]-a_positive[p]))*(log(2)/halflifeTimeDifference));
            sumCouplePositive = sumCouplePositive + couple;
            couple = 0;
          }
        }
        if (a_negative.size() >= i_negative.size() && i_negative.size() != 0) {
          for (int p = 0; p < i_negative.size(); p++){
            couple = 1 * exp(-(std::abs(i_negative[p]-a_negative[p]))*(log(2)/halflifeTimeDifference));
            sumCoupleNegative = sumCoupleNegative + couple;
            couple = 0;
          }
        }
        if (i_positive.size() > a_positive.size() && a_positive.size() != 0) {
          for (int p = 0; p < a_positive.size(); p++){
            couple = 1 * exp(-(std::abs(i_positive[p]-a_positive[p]))*(log(2)/halflifeTimeDifference));
            sumCouplePositive = sumCouplePositive + couple;
            couple = 0;
          }
        }
        if (i_negative.size() > a_negative.size() && a_negative.size() != 0) {
          for (int p = 0; p < a_negative.size(); p++){
            couple = 1 * exp(-(std::abs(i_negative[p]-a_negative[p]))*(log(2)/halflifeTimeDifference));
            sumCoupleNegative = sumCoupleNegative + couple;
            couple = 0;
          }
        }
        //// how many actors used concept v[k] in same manner as a? // how many concepts are used by both in same manner?
        sumConcept = sumConcept + sumCoupleNegative + sumCouplePositive;
      }else{ // if matchNomatchSim = "nomatch"
             // if there acctually is an intersection in xw
             if (xw.size() != 0 ) {
               
               // for each entry in xw:
                 for (int m = 0; m < xw.size(); m++) {
                   i_sendertarget.clear();
                   a_sendertarget.clear();
                   
                   // loop back over all events until i
                   for (int n = 0; n < i; n++){
                     if ( senderTargetSim == "sender" ){
                       if (sender[n] == v[k] && target[n] == xw[m] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                         i_sendertarget.push_back(time[n]);
                       }
                       if (sender[n] == sender[i] && target[n] == xw[m]  && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                         a_sendertarget.push_back(time[n]);
                       }
                     }else{
                       if (target[n] == v[k] && sender[n] == xw[m] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                         i_sendertarget.push_back(time[n]);
                       }
                       if (target[n] == target[i] && sender[n] == xw[m] && time[n] != time[i] && eventAttributeVar[n] == eventAttribute){
                         a_sendertarget.push_back(time[n]);
                       }
                     }//closes senderTargetSim == target
                   }//closes n-loop
                 }//closes m-loop
             }//closes xw.size() != 0
             
             sumCoupleNoMatch = 0;
             
             //how large are the respextive vectors with the matches in them?
             if (a_sendertarget.size() >= i_sendertarget.size() && i_sendertarget.size() != 0) {
               for (int p = 0; p < i_negative.size(); p++){
                 couple = 1 * exp(-(std::abs(i_sendertarget[p]-a_sendertarget[p]))*(log(2)/halflifeTimeDifference));
                 sumCoupleNoMatch = sumCoupleNoMatch + couple;
                 couple = 0;
               }
             }
             if (i_sendertarget.size() > a_sendertarget.size() && a_sendertarget.size() != 0) {
               for (int p = 0; p < a_positive.size(); p++){
                 couple = 1 * exp(-(std::abs(i_sendertarget[p]-a_sendertarget[p]))*(log(2)/halflifeTimeDifference));
                 sumCoupleNoMatch = sumCoupleNoMatch + couple;
                 couple = 0;
               }
             }
             //add them together
             sumConcept = sumConcept + sumCoupleNoMatch;	
      }//closes matchNomatchSim == "nomatch"
      
      // how many actors used concept v[k]? // how many concepts are used by both?
      if ( xw.size() != 0 ){
        // find time for the time-discount in the simple-similarity equation
        if ( senderTargetSim == "sender"){
          // find time, when actor k used concept $p$ last
          for (int q = i-1; q >= 0; q--){
            if ( matchNomatchSim == "nomatch"){
              if (sender[q] == v[k] && target[q] == target[i] && time[q] != time[i] &&
                    eventAttributeVar[q] == eventAttribute){
                timePLast = time[q];
                break;
              }
            }else if (matchNomatchSim == "match"){
              if (sender[q] == v[k] && target[q] == target[i] && time[q] != time[i] &&
                    eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == eventTypeVar[i]){
                timePLast = time[q];
                break;
              }
            }
          }//closes q-loop
        }else if ( senderTargetSim == "target"){
          // find time, when target k was last used by $a$
            for (int q = i-1; q >= 0; q--){
              if (matchNomatchSim == "nomatch"){
                if (target[q] == v[k] && sender[q] == sender[i] && time[q] != time[i] &&
                      eventAttributeVar[q] == eventAttribute){
                  timePLast = time[q];
                  break;
                }
              } else if (matchNomatchSim == "match"){
                if (target[q] == v[k] && sender[q] == sender[i] && time[q] != time[i] &&
                      eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == eventTypeVar[i]){
                  timePLast = time[q];
                  break;
                }
              }
            }//closes q-loop
        }//closes senderTargetSim == "target"
        
        // calculate weight of each actor/target
        weightSim = sumConcept * exp(-(time[i]-timePLast)*xlog) * xlog;
        totalSim = totalSim + weightSim;
        if (weightSim != 0) {
          counter++;   
        }
      }//closes if xw.size != 0
    }//closes k-loop
    
    if (counter == 0) {
      result[i] = 0;
    }else{
      result[i] = (totalSim/counter);
    }
  }//i-loop close
  return Rcpp::wrap(result);
}

