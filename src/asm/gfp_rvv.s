  .text
  .balign 4
  .global gfp_neg
  .global gfp_add
  .global gfp_sub
  .global gfp_mul

# void gfp_carry(uint64_t *p2, uint64_t *c, uint64_t head);
# a0 => *p2
# a1 => *c
# a2 => head
# gfp_carry:
#   li t0, 4                          # Load 4 uint64_t integer
#   vsetvli t1, t0, e64, m2, ta, ma   # fixed length vectors of elements (qemu with 128 bit register)
#   vle64.v v2, (a0)                  # Load p2
#   vle64.v v4, (a1)                  # Load c
#   vsbc.vvm v2, v4, v2, v0           # [code] bi = ci - pi - carry
#   sub a3, x0, a2                    # [code] a3 = 0 - head - carry
#   bnez a3, 1f                       # [code] if (a3 == 0) { *c = b }
#   vse64.v v2, (a1)                  # [code] *c = b
# 1:
#   ret

# void gfp_neg(uint64_t *p2, uint64_t *c, uint64_t *a);
# a0 => *p2
# a1 => *c
# a2 => *a
gfp_neg:
  li t0, 4
  vsetvli t1, t0, e64, m2, ta, ma  # fixed length vectors of elements (qemu with 128 bit register)
  vle64.v v2 (a0)                  # Load p2
  vle64.v v4 (a2)                  # Load a
  vsbc.vmm v4, v2, v4, v0          # [code] ci(v4) = pi(v2) - ai(v4) - carry
  # gfp_carry()
  vsbc.vmm v2, v4, v2, v0          # [code] bi(v2) = ci(v4) - pi(v2) - carry
  li a3, 0                         # [code] a3 = 0 - head - carry
  bnez a3, 1f                      # [code] if (a3 == 0) { *c = b }
  vse64.v v2, (a1)
1:
  ret

# void gfp_add(uint64_t *p2, uint64_t *c, uint64_t *a, uint64_t *b);
# a0 => *p2
# a1 => *c
# a2 => *a
# a3 => *b
gfp_add:
  li t0, 4
  vsetvli t1, t0, e64, m2, ta, ma  # fixed length vectors of elements (qemu with 128 bit register)
  vle64.v v2 (a2)                  # Load a
  vle64.v v4 (a3)                  # Load b
  vadc.vmm v4, v2, v4, v0          # [code] ci = ai + bi + carry
  # gfp_carry()
  vle64.v v2 (a0)                  # Load p2
  vsbc.vmm v2, v4, v2, v0          # [code] bi(v2) = ci(v4) - pi(v2) - carry
  li a3, 0                         # [code] a3 = 0 - head - carry
  bnez a3, 1f                      # [code] if (a3 == 0) { *c = b }
  vse64.v v2, (a1)
1:
  ret

# void gfp_sub(uint64_t *p2, uint64_t *c, uint64_t *a, uint64_t *b);
# a0 => *p2
# a1 => *c
# a2 => *a
# a3 => *b
gfp_sub:
  li t0, 4
  vsetvli t1, t0, e64, m2, ta, ma  # fixed length vectors of elements (qemu with 128 bit register)
  vle64.v v2 (a0)                  # Load p2
  vle64.v v4 (a2)                  # Load a
  vle64.v v6 (a3)                  # Load b
  vsbc.vmm v6, v2, v6, v0          # [code] ti(v6) = pi(v2) - bi(v6) - carry
  vadc.vmm v4, v4, v6, v0          # [code] ci(v4) = ai(v4) + ti(v6) + carry
  # gfp_carry()
  vsbc.vmm v2, v4, v2, v0          # [code] bi(v2) = ci(v4) - pi(v2) - carry
  li a3, 0                         # [code] a3 = 0 - head - carry
  bnez a3, 1f                      # [code] if (a3 == 0) { *c = b }
  vse64.v v2, (a1)
  ret

# void gfp_mul(uint64_t *p2, uint64_t *c, uint64_t *a, uint64_t *b);
# a0 => *p2
# a1 => *c
# a2 => *a
# a3 => *b
gfp_mul:
  ret
