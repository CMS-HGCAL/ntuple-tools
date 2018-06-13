//
//  Helpers.hpp
//  xHGCalClustering
//
//  Created by Jeremi Niedziela on 13/06/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Helpers_h
#define Helpers_h

bool pointWithinCircle(double px,double py,double x,double y,double r){
  return pow(r,2) >= pow(px-x,2) + pow(py-y,2);
}


double duration(std::chrono::time_point<std::chrono::system_clock> t0,
                std::chrono::time_point<std::chrono::system_clock> t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

std::chrono::time_point<std::chrono::system_clock> now()
{
  return std::chrono::system_clock::now();
}

#endif /* Helpers_h */
