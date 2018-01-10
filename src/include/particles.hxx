
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"

// ======================================================================
// mparticles_base

template<typename S>
struct mparticles_base
{
  using sub_t = S;

  explicit mparticles_base(psc_mparticles *mprts) : mprts_(mprts) { }

  void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
  {
    psc_mparticles_put_as(mprts_, mprts_base, flags);
  }

  unsigned int n_patches() { return mprts_->nr_patches; }
  
  psc_mparticles *mprts() { return mprts_; }
  
  sub_t* sub() { return mrc_to_subobj(mprts(), sub_t); }

private:
  psc_mparticles *mprts_;
};

// ======================================================================
// psc_particle_buf

template<typename P>
struct psc_particle_buf
{
  using particle_t = P;

  psc_particle_buf()
    : m_data(), m_size(), m_capacity()
  {
  }

  ~psc_particle_buf()
  {
    free(m_data);
  }

  psc_particle_buf(const psc_particle_buf&) = delete;
  
  void resize(unsigned int new_size)
  {
    assert(new_size <= m_capacity);
    m_size = new_size;
  }

  void reserve(unsigned int new_capacity)
  {
    if (new_capacity <= m_capacity)
      return;

    new_capacity = std::max(new_capacity, m_capacity * 2);
    
    m_data = (particle_t *) realloc(m_data, new_capacity * sizeof(*m_data));
    m_capacity = new_capacity;
  }

  void push_back(const particle_t& prt)
  {
    unsigned int n = m_size;
    if (n >= m_capacity) {
      reserve(n + 1);
    }
    m_data[n++] = prt;
    m_size = n;
  }

  particle_t& operator[](int n)
  {
    return m_data[n];
  }

  particle_t *m_data;
  unsigned int m_size;
  unsigned int m_capacity;

  unsigned int size() const { return m_size; }
  unsigned int capacity() const { return m_capacity; }
};

// ======================================================================
// mparticles

template<typename S>
struct mparticles : mparticles_base<S>
{
  using Base = mparticles_base<S>;
  using particles_t = typename Base::sub_t::particles_t;
  
  mparticles(psc_mparticles *mprts) : mparticles_base<S>(mprts) { }

  particles_t& operator[](int p)
  {
    return this->sub()->patch[p];
  }

};

#endif

