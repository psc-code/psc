
#pragma once

// ======================================================================
// PscCollision

template<typename S>
struct PscCollision
{
  using sub_t = S;
  
  explicit PscCollision(psc_collision *collision)
    : collision_(collision)
  {}
  
  sub_t* sub() { return mrc_to_subobj(collision_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_collision *collision_;
};

