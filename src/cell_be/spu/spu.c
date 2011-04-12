int spu_main(unsigned long long spe_id, unsigned long long spu_comm_ea,
	     unsigned long long ea);

int
main(unsigned long long spe_id, unsigned long long spu_comm_ea,
     unsigned long long env)
{
  return spu_main(spe_id, spu_comm_ea, env);
}
