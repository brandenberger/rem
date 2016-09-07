#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;


//####################################################################
// [[Rcpp::export]]
NumericVector absoluteDiffAverageWeightEventAttributeCpp(
    std::vector<std::string> sender,
    std::vector<std::string> target,
    NumericVector time,
    NumericVector weightvar,
    std::vector<std::string> eventattributevar,
    std::string eventattribute,
    double xlog) {
    
    NumericVector result(time.size());
    int count = 0;
    
    // for-loop: for each event, do:
      for ( int i = 0; i < time.size(); i++){
        
        double totaldiff = 0;
        count = 0;
        
        for ( int w = 0; w < i; w++ ) {
          if (eventattributevar[w] == eventattribute && sender[i] != sender[w] && target[i] == target[w]) {
            if ( time[i] != time[w] ) {
              count = count + 1;
              totaldiff = totaldiff + std::abs(weightvar[i] - weightvar[w]);
            } 
          }
        } // closes w-loop
        if (count == 0) {
          result[i] = 0;
        }else{
          result[i] = totaldiff/count;
        }
      } // closes i-loop
    return Rcpp::wrap(result);
}


//####################################################################
// [[Rcpp::export]]
NumericVector triadOldCpp(
  std::vector<std::string> sender, 
  std::vector<std::string> target,
  NumericVector time, 
  NumericVector weightvar, 
  std::vector<std::string> typevar, 
  std::string typeA, 
  std::string typeB,
  std::vector<std::string> attributevarAB, 
  std::string attrAB,
  std::vector<std::string> attributevarAI,
  std::string attrAI,
  std::vector<std::string> attributevarBI,
  std::string attrBI,
  double xlog) {
    
    NumericVector result(sender.size());
    double weighta;
    double weightb;
    double totalweighta;
    double totalweightb;
    double weightab;
    double totalweight;
    std::vector<std::string> x;
    std::vector<std::string> y;
    std::vector<std::string> v;
    
    for ( int i = 0; i < sender.size(); i++){
      
      //filter out events i that have attribute attrAB
      if ( attributevarAB[i] == attrAB ){
        
        weightab = 0;
        totalweight = 0;
        // clear strings
        x.clear();
        y.clear();
        v.clear();
        
        // get list of partners s and t have interacted with
        for ( int w = 0; w < i-1; w++ ) {       
          // with whom has sender interacted in the past?
          if ( sender[w] == sender[i] && target[w] != target[i] && 
          time[w] != time[i] && typevar[w] == typeA && attributevarAI[w] == attrAI ){
            x.push_back(target[w]);
          }
          if ( target[w] == sender[i] && sender[w] != target[i] &&
          time[w] != time[i] && typevar[w] == typeA && attributevarAI[w] == attrAI ){
            x.push_back(sender[w]);
          }
          // with whom has target interacted in the past?
          if ( sender[w] == target[i] && target[w] != sender[i] && 
          time[w] != time[i] &&  typevar[w] == typeB && attributevarBI[w] == attrBI ){
            y.push_back(target[w]);
          }
          if ( target[w] == target[i] && sender[w] != sender[i] &&
          time[w] != time[i] && typevar[w] == typeB && attributevarBI[w] == attrBI ){
            y.push_back(sender[w]);
          }
        } // close w-loop
        
        //order x and y - and create v: the array that holds all the actors a and 
        //b have interacted with in the past
        std::sort(x.begin(),x.end());
        std::sort(y.begin(),y.end());
        // find intersection between y and x string entries (http://www.cplusplus.com/reference/algorithm/set_intersection/)
        std::set_intersection (x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(v) );                                              
        // remove duplicates from v
        sort( v.begin(), v.end() );
        v.erase( unique( v.begin(), v.end() ), v.end() );
        
        // for each entry in v
        for (size_t j = 0; j < v.size(); j++) {
          totalweighta = 0;
          totalweightb = 0;
          for ( int z = 0; z < i-1; z++ ) {   
            weighta = 0;
            weightb = 0;
            //caluculate weighta
            if ( ( (sender[z] == sender[i] && target[z] == v[j]) || 
            (target[z] == sender[i] && sender[z] == v[j]) ) && typevar[z] == typeA &&
            attributevarAI[z] == attrAI && time[z] != time[i] ){
              weighta = std::abs(weightvar[z]) * exp( - ( time[i] - time[z] ) * xlog) * xlog;
              totalweighta = totalweighta + weighta;   
            }
            //calculate weightb
            if ( ( (sender[z] == target[i] && target[z] == v[j]) || 
            (target[z] == target[i] && sender[z] == v[j]) ) && typevar[z] == typeB &&
            attributevarBI[z] == attrBI && time[z] != time[i] ){
              weightb = std::abs(weightvar[z]) * exp( - ( time[i] - time[z] ) * xlog) * xlog;
              totalweightb = totalweightb + weightb;   
            }
          } //closes z-loop
          
          // multiply totalweighta times totalweightb
          weightab = totalweighta * totalweightb;
          totalweight = totalweight + weightab;   
        } //closes j-loop  
        //take the squared rood of totalweight
        result[i] = sqrt(totalweight);
      }//closes if-attrvar == attrAB
    } // closes i-loop
    return Rcpp::wrap(result);
  }













