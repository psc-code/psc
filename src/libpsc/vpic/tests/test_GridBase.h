
#ifndef TEST_GRID_BASE_H
#define TEST_GRID_BASE_H

template<typename Grid>
Grid* test_GridBase_create()
{
  int gdims[3] = { 16, 32, 1 };
  double xl[3] = { -1., -2., -4. };
  double xh[3] = {  1. , 2. , 4. };
  int np[3] = { 1, 1, 1 };
  double dt = .05;
  double cvac = 1.;
  double eps0 = 1.;

  double dx[3];
  for (int d = 0; d < 3; d++) {
    dx[d] = (xh[d] - xl[d]) / gdims[d];
  }

  Grid *grid = Grid::create();
  grid->setup(dx, dt, cvac, eps0);
  grid->partition_periodic_box(xl, xh, gdims, np);

  return grid;
}

template<typename Grid>
void test_GridBase_methods(Grid* g)
{
  (void) g->mp_send_buffer(1);
  (void) g->mp_recv_buffer(1);
}

template<typename Grid>
void test_GridBase_destroy(Grid* g)
{
}

template<typename Grid>
void test_GridBase()
{
  Grid *g = test_GridBase_create<Grid>();
  test_GridBase_methods(g);
  test_GridBase_destroy(g);
}


#endif
