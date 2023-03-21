
#ifdef INLIB_MEM
#include <inlib/mem>
#endif //INLIB_MEM

#include <inlib/lina/umat>

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//#include <inlib/lina/matTs>
namespace inlib {
////////////////////////////////////////////////
/// D=3 ////////////////////////////////////////
////////////////////////////////////////////////
template <class MAT>
inline void set_R1(MAT& a_m) {
  matrix_set<MAT>(a_m, 0, 0, 0,
                       0, 0, 1,
                       0,-1, 0);
}
template <class MAT>
inline void set_R2(MAT& a_m) {
  matrix_set<MAT>(a_m, 0, 0,-1,
                       0, 0, 0,
                       1, 0, 0);
}
template <class MAT>
inline void set_R3(MAT& a_m) {
  matrix_set<MAT>(a_m, 0, 1, 0,
                      -1, 0, 0,
                       0, 0, 0);
}

template <class VECTOR>
inline void set_rots(VECTOR& a_rots) {
  a_rots.resize(3);
  set_R1(a_rots[0]);
  set_R2(a_rots[1]);
  set_R3(a_rots[2]);
}
////////////////////////////////////////////////
/// D=4 ////////////////////////////////////////
////////////////////////////////////////////////
template <class MAT>
inline void set_eta(MAT& a_m) {
  a_m.set_zero();
  a_m.set_value(0,0, 1);
  a_m.set_value(1,1,-1);
  a_m.set_value(2,2,-1);
  a_m.set_value(3,3,-1);
}

template <class MAT>
inline void set_J(MAT& a_m) {
  matrix_set(a_m,    0, 0, 0, 1,
                     0, 0,-1, 0,
                     0, 1, 0, 0,
                    -1, 0, 0, 0);
}
template <class MAT>
inline void set_E1(MAT& a_m) {
  matrix_set(a_m,        0, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 1,
                         0, 0,-1, 0);
}

template <class MAT>
inline void set_E2(MAT& a_m) {
  matrix_set(a_m,        0, 0, 0, 0,
                         0, 0, 0,-1,
                         0, 0, 0, 0,
                         0, 1, 0, 0);
}

template <class MAT>
inline void set_E3(MAT& a_m) {
  matrix_set(a_m,        0, 0, 0, 0,
                         0, 0, 1, 0,
                         0,-1, 0, 0,
                         0, 0, 0, 0);
}

template <class MAT>
inline void set_E4(MAT& a_m) {
  matrix_set(a_m,        0, 0, 0, 1,
                         0, 0, 0, 0,
                         0, 0, 0, 0,
                        -1, 0, 0, 0);
}

template <class MAT>
inline void set_E5(MAT& a_m) {
  matrix_set(a_m,        0, 0, 1, 0,
                         0, 0, 0, 0,
                        -1, 0, 0, 0,
                         0, 0, 0, 0);
}

template <class MAT>
inline void set_E6(MAT& a_m) {
  matrix_set(a_m,        0, 1, 0, 0,
                        -1, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 0);
}
template <class VECTOR>
inline void set_Es(VECTOR& a_Es) {
  //  | 0  0  0  0|  | 0  0  0  0|  | 0  0  0  0|  | 0  1  0  0|  | 0  0  1  0|  | 0  0  0  1|
  //  | 0  0  0  0|  | 0  0  0 -1|  | 0  0  1  0|  |-1  0  0  0|  | 0  0  0  0|  | 0  0  0  0|
  //  | 0  0  0  1|  | 0  0  0  0|  | 0 -1  0  0|  | 0  0  0  0|  |-1  0  0  0|  | 0  0  0  0|
  //  | 0  0 -1  0|  | 0  1  0  0|  | 0  0  0  0|  | 0  0  0  0|  | 0  0  0  0|  |-1  0  0  0|
  a_Es.resize(6);
  set_E1(a_Es[0]);
  set_E2(a_Es[1]);
  set_E3(a_Es[2]);
  //WARNING E6,5,4 and note E4,5,6
  set_E6(a_Es[3]);
  set_E5(a_Es[4]);
  set_E4(a_Es[5]);
}
template <class VECTOR>
inline void set_Es_eta(VECTOR& a_vec) {
  typedef typename VECTOR::value_type MAT;
  a_vec.clear();
  VECTOR Es;
  set_Es(Es);
  MAT eta;
  set_eta(eta);
 {for(size_t a=0;a<6;a++) a_vec.push_back(Es[a]*eta);}
}

////////////////////////////////////
/// i*D_majoranas : ////////////////
////////////////////////////////////
template <class MAT>
inline void set_M0(MAT& a_m) {
  set_J(a_m);
}

template <class MAT>
inline void set_M1(MAT& a_m) {
  matrix_set(a_m,   -1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0,-1, 0,
                     0, 0, 0, 1);
}

template <class MAT>
inline void set_M2(MAT& a_m) {
  matrix_set(a_m,    0, 0, 0,-1,
                     0, 0, 1, 0,
                     0, 1, 0, 0,
                    -1, 0, 0, 0);
}

template <class MAT>
inline void set_M3(MAT& a_m) {
  matrix_set(a_m,    0, 1, 0, 0,
                     1, 0, 0, 0,
                     0, 0, 0, 1,
                     0, 0, 1, 0);
}

template <class VECTOR>
inline void set_majoranas_up(VECTOR& a_majs) {
  a_majs.resize(4);
  set_M0(a_majs[0]);
  set_M1(a_majs[1]);
  set_M2(a_majs[2]);
  set_M3(a_majs[3]);
  typedef typename VECTOR::value_type MAT; //mat4
  typedef typename MAT::elem_t T; //std::complex
  T minus_zi(0,-1);
  for(size_t mu=0;mu<4;mu++) a_majs[mu] *= minus_zi;
}

template <class VECTOR>
inline void set_majoranas_down(VECTOR& a_majs) {
  set_majoranas_up(a_majs);
  for(size_t j=0;j<3;j++) a_majs[j+1].multiply(-1);
}

template <class MAT>
inline void set_maj_0_S(MAT& a_m) {
  matrix_set(a_m,    1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);
}
template <class MAT>
inline void set_maj_1_S(MAT& a_m) {
  matrix_set(a_m,    0, 0, 0,-1,
                     0, 0,-1, 0,
                     0,-1, 0, 0,
                    -1, 0, 0, 0);
}
template <class MAT>
inline void set_maj_2_S(MAT& a_m) {
  matrix_set(a_m,    1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0,-1, 0,
                     0, 0, 0,-1);
}
template <class MAT>
inline void set_maj_3_S(MAT& a_m) {
  matrix_set(a_m,    0, 0,-1, 0,
                     0, 0, 0, 1,
                    -1, 0, 0, 0,
                     0, 1, 0, 0);
}

template <class MAT>
inline void set_maj_D1_S(MAT& a_m) {
  matrix_set(a_m,    0,-1, 0, 0,
                    -1, 0, 0, 0,
                     0, 0, 0, 1,
                     0, 0, 1, 0);
  a_m.multiply(0.5);
}
template <class MAT>
inline void set_maj_D2_S(MAT& a_m) {
  matrix_set(a_m,    0, 0,-1, 0,
                     0, 0, 0,-1,
                    -1, 0, 0, 0,
                     0,-1, 0, 0);
  a_m.multiply(0.5);
}
template <class MAT>
inline void set_maj_D3_S(MAT& a_m) {
  matrix_set(a_m,   -1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0,-1);
  a_m.multiply(0.5);
}
template <class MAT>
inline void set_maj_D4_S(MAT& a_m) {
  matrix_set(a_m,    1, 0, 0, 0,
                     0,-1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0,-1);
  a_m.multiply(0.5);
}
template <class MAT>
inline void set_maj_D5_S(MAT& a_m) {
  matrix_set(a_m,    0, 0, 0, 1,
                     0, 0,-1, 0,
                     0,-1, 0, 0,
                     1, 0, 0, 0);
  a_m.multiply(0.5);
}
template <class MAT>
inline void set_maj_D6_S(MAT& a_m) {
  matrix_set(a_m,    0,-1, 0, 0,
                    -1, 0, 0, 0,
                     0, 0, 0,-1,
                     0, 0,-1, 0);
  a_m.multiply(0.5);
}
////////////////////////////////////////////////
/// D=5 ////////////////////////////////////////
////////////////////////////////////////////////
template <class MAT>
inline void set_eta_131(MAT& a_m) {
  matrix_set(a_m,    1, 0, 0, 0, 0,
                     0,-1, 0, 0, 0,
                     0, 0,-1, 0, 0,
                     0, 0, 0,-1, 0,
                     0, 0, 0, 0, 1);
}

template <class VECTOR>
inline void set_Es_5(VECTOR& a_Es) { //10 matrices.
  a_Es.resize(10);

  matrix_set(a_Es[0],    0, 0, 0, 0, 1,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                        -1, 0, 0, 0, 0);
  matrix_set(a_Es[1],    0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0,-1, 0, 0, 0);
  matrix_set(a_Es[2],    0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1,
                         0, 0, 0, 0, 0,
                         0, 0,-1, 0, 0);
  matrix_set(a_Es[3],    0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1,
                         0, 0, 0,-1, 0);

  matrix_set(a_Es[4],    0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 1, 0,
                         0, 0,-1, 0, 0,
                         0, 0, 0, 0, 0);
  matrix_set(a_Es[5],    0, 0, 0, 0, 0,
                         0, 0, 0,-1, 0,
                         0, 0, 0, 0, 0,
                         0, 1, 0, 0, 0,
                         0, 0, 0, 0, 0);
  matrix_set(a_Es[6],    0, 0, 0, 0, 0,
                         0, 0, 1, 0, 0,
                         0,-1, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0);

  matrix_set(a_Es[7],    0, 1, 0, 0, 0,
                        -1, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0);
  matrix_set(a_Es[8],    0, 0, 1, 0, 0,
                         0, 0, 0, 0, 0,
                        -1, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0);
  matrix_set(a_Es[9],    0, 0, 0, 1, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                        -1, 0, 0, 0, 0,
                         0, 0, 0, 0, 0);
}
////////////////////////////////////////////////
/// D=6 ////////////////////////////////////////
////////////////////////////////////////////////
template <class MAT>
inline void set_G6(MAT& a_m) {
  //    | I3   0|
  //    | 0  -I3|
  a_m.set_zero();
  a_m.set_value(0,0,1);
  a_m.set_value(1,1,1);
  a_m.set_value(2,2,1);
  a_m.set_value(3,3,-1);
  a_m.set_value(4,4,-1);
  a_m.set_value(5,5,-1);
}

template <class MAT>
inline void set_H6(MAT& a_m) {
  //    | 0  I3|
  //    | I3  0|
  a_m.set_zero();
  a_m.set_value(0,3,1);
  a_m.set_value(1,4,1);
  a_m.set_value(2,5,1);
  a_m.set_value(3,0,1);
  a_m.set_value(4,1,1);
  a_m.set_value(5,2,1);
}

////////////////////////////////////////////////
/// D=10 ///////////////////////////////////////
////////////////////////////////////////////////
template <class MAT>
inline void set_eta_G6(MAT& a_m) {
  //    | eta   0|
  //    | 0    G6|

  typedef typename MAT::elem_t T;
  T minus_one = T(-1);

  a_m.set_zero();
  a_m.set_value(0,0, 1);
  a_m.set_value(1,1, minus_one);
  a_m.set_value(2,2, minus_one);
  a_m.set_value(3,3, minus_one);
  a_m.set_value(4,4, 1);
  a_m.set_value(5,5, 1);
  a_m.set_value(6,6, 1);
  a_m.set_value(7,7, minus_one);
  a_m.set_value(8,8, minus_one);
  a_m.set_value(9,9, minus_one);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
template <size_t NUM,class VECTOR>
inline void set_from_csts_down(VECTOR& a_Cs,
                          const typename VECTOR::value_type::elem_t a_csts[NUM][NUM][NUM],
                          typename VECTOR::value_type::elem_t a_factor = 1) {
  typedef typename VECTOR::value_type REP;
  typedef typename REP::elem_t T;
  a_Cs.resize(NUM);
 {for(size_t a=0;a<NUM;a++) {
    REP& _mat = a_Cs[a];
    for(size_t j=0;j<NUM;j++) {
      for(size_t k=0;k<NUM;k++) {
        _mat.set_value(j,k,a_factor*a_csts[a][k][j]); // with the "down" convention.
      }
  }}}
}

template <class VECTOR,class VECA>
inline bool set_rot_boost_A_rep(const VECTOR& a_diracs,const VECA& a_As,double a_factor,VECTOR& a_reps) {
  // a_diracs : for exa the four dirac matrices.
  // a_As : for example Es.
  // a_reps will have a_As.size() matrices : a_rep[a] = a_factor*[a_diracs[mu],a_diracs[nu]]*a_As[a](mu,nu)
  typedef typename VECTOR::value_type MAT;
  typedef typename MAT::elem_t T;
  size_t _sz = a_As.size();
  if(!_sz) return false;
  a_reps.resize(_sz);
  size_t _dim = a_As[0].dimension();
  MAT tmp,cm;
  for(size_t a=0;a<_sz;a++) {
    MAT& res = a_reps[a];
    res.set_zero();
    for(size_t mu=0;mu<_dim;mu++) {
    for(size_t nu=0;nu<_dim;nu++) {
      commutator(a_diracs[mu],a_diracs[nu],tmp,cm);
      cm.multiply(a_As[a].value(mu,nu));
      res += cm;
    }}
    res.multiply(T(a_factor));
  }
  return true;
}

template <class VECTOR,class MATRIX_NUM,size_t NUM,class PREC>
inline bool set_commutator_constants_prec(const VECTOR& a_Ms,
  typename VECTOR::value_type::elem_t a_csts[NUM][NUM][NUM],
  const PREC& a_prec,PREC(*a_fabs)(const typename VECTOR::value_type::elem_t&)
) {
  // If exists: [Ms(a),Ms(b)] = a_csts[a][b][c]*Ms(c)
  // If so, we have:
  //   Trace[Ms(d)*[Ms(a),Ms(b)]] = a_csts[a][b][c]*Trace[Ms(d)*Ms(c)]
  // If inv of Trace[Ms(d)*Ms(c)] (d=row,c=col) exists:
  //   a_csts[a][b][e] = Trace[Ms(e)*Ms(d)]_inv * Trace[Ms(d)*[Ms(a),Ms(b)]]

  typedef typename VECTOR::value_type MATRIX;
  typedef typename MATRIX::elem_t T;

  MATRIX_NUM tr_Ms_Ms;
  MATRIX tmp;
 {for(size_t a=0;a<NUM;a++) {
  for(size_t b=0;b<NUM;b++) {
    tmp = a_Ms[a]*a_Ms[b];
    tr_Ms_Ms.set_value(a,b,tmp.trace());
  }}}
  //mat_dump(std::cout,"tr_Ms_Ms",tr_Ms_Ms);
  MATRIX_NUM tr_Ms_Ms_inv;
  if(tr_Ms_Ms.is_diagonal_prec(a_prec,a_fabs)) {
    for(size_t a=0;a<NUM;a++) {
      T value = tr_Ms_Ms.value(a,a);
      if(value==T(0)) return false;
      tr_Ms_Ms_inv.set_value(a,a,T(1)/value);
    }
  } else {
    //mat_dump(std::cout,"tr_Ms_Ms not diagonal",tr_Ms_Ms);
    //the below is too long with NUM=16. Have a faster invert_sym() method ?
    if(T(tr_Ms_Ms.invert(tr_Ms_Ms_inv))==T(0)) return false;
  }
  //mat_dump(std::cout,"tr_Ms_Ms_inv",tr_Ms_Ms_inv);

  //   a_csts[a][b][e] = Trace[Ms(e)*Ms(d)]_inv * Trace[Ms(d)*[Ms(a),Ms(b)]]
 {for(size_t a=0;a<NUM;a++) {
  for(size_t b=0;b<NUM;b++) {
  for(size_t e=0;e<NUM;e++) {
    T right = T(0);
    for(size_t d=0;d<NUM;d++) {
      tmp = a_Ms[d]*(a_Ms[a]*a_Ms[b]-a_Ms[b]*a_Ms[a]);
      T tmp_trace = tmp.trace();
      right += tr_Ms_Ms_inv(e,d)*tmp.trace();
    }
    a_csts[a][b][e] = right;
  }}}}

  return true;
}

}


//#include <inlib/lina/Tgroup>
namespace inlib {
template <class T>
inline void set_Es_eta_group_constants(T a_a[6][6][6]) {
  // Si M[a].set_value(j,k) = csts[a][j][k]  // with the "up" convention.
  // M 6 6x6 matrices :
  //  M123 = |   e123     0 |
  //         |      0  e123 |
  //  M456 = |      0  e123 |
  //         |  -e123     0 |
  //
  // adjs :
  //   |-e123      0|    |    0 -e123|
  //   |   0   -e123|    | e123     0|
  // with : e123  e123=e3

  // e1 | 0  0  0 |  e2 | 0  0 -1 | e3 | 0  1  0 |
  //    | 0  0  1 |     | 0  0  0 |    |-1  0  0 |
  //    | 0 -1  0 |     | 1  0  0 |    | 0  0  0 |

  //E1,E2,E3,E6,E5,E4
  // 0  1  2  3  4  5

  // We must have : [Es(a)*eta,Es(b)*eta] = a_a(a,b,c) * Es(c)*eta

  for(unsigned int i=0;i<6;i++) {
    for(unsigned int j=0;j<6;j++) {
      for(unsigned int k=0;k<6;k++) {
        a_a[i][j][k] = T();
      }
    }
  }

  a_a[0][1][2] = 1;
  a_a[0][2][1] = -1;
  a_a[0][4][5] = 1;
  a_a[0][5][4] = -1;
  a_a[1][0][2] = -1;
  a_a[1][2][0] = 1;
  a_a[1][3][5] = -1;
  a_a[1][5][3] = 1;
  a_a[2][0][1] = 1;
  a_a[2][1][0] = -1;
  a_a[2][3][4] = 1;
  a_a[2][4][3] = -1;
  a_a[3][1][5] = 1;
  a_a[3][2][4] = -1;
  a_a[3][4][2] = -1;
  a_a[3][5][1] = 1;
  a_a[4][0][5] = -1;
  a_a[4][2][3] = 1;
  a_a[4][3][2] = 1;
  a_a[4][5][0] = -1;
  a_a[5][0][4] = 1;
  a_a[5][1][3] = -1;
  a_a[5][3][1] = -1;
  a_a[5][4][0] = 1;
}
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#include <inlib/mathd>
#include <inlib/mathz>
#include <inlib/random>

#include <cstdio>
#include <cstdlib>
#include <iostream>

int main(int argc,char** argv) {
#ifdef INLIB_MEM
  inlib::mem::set_check_by_class(true);{
#endif //INLIB_MEM

  bool verbose = (argc==2 && !::strcmp(argv[1],"-verbose"))?true:false;

  double depsilon = 1e-9;
//inlib::random::gauss rg;
  inlib::random::gauss rg(1,2);
  
  typedef std::complex<double> z;
  z zi(0,1);

  //////////////////////////////////////////////////////////
  /// D=3 material: ////////////////////////////////////////
  //////////////////////////////////////////////////////////
  typedef inlib::mat<double,3> m3d;
  std::vector<m3d> rs;inlib::set_rots(rs);

  //////////////////////////////////////////////////////////
  /// D=4 material: ////////////////////////////////////////
  //////////////////////////////////////////////////////////
  typedef inlib::mat<double,4> m4d;
  
  m4d eta;inlib::set_eta(eta);  // the eta metric matrix.
  std::vector<m4d> Es;inlib::set_Es(Es); //The A[a] in the SMAD paper.
  if(verbose) mat_dump(std::cout,"eta:",eta);
  
  std::vector<m4d> Es_eta;inlib::set_Es_eta(Es_eta);
  std::vector<m4d> minus_Es_eta(Es_eta);
  inlib::multiply(minus_Es_eta,-1.0);

  std::vector<m4d> tilde_Ss(4);   //4 S_tild matrices of the SMAD paper.
  inlib::set_maj_0_S(tilde_Ss[0]);
  inlib::set_maj_1_S(tilde_Ss[1]);
  inlib::set_maj_2_S(tilde_Ss[2]);
  inlib::set_maj_3_S(tilde_Ss[3]);

  std::vector<m4d> DSs(6); //6 S matrices of the SMAD paper.
  inlib::set_maj_D1_S(DSs[0]);
  inlib::set_maj_D2_S(DSs[1]);
  inlib::set_maj_D3_S(DSs[2]);
  inlib::set_maj_D4_S(DSs[3]);
  inlib::set_maj_D5_S(DSs[4]);
  inlib::set_maj_D6_S(DSs[5]);

  std::vector<m4d> SIGMAs; //10 SIGMAs marices of the SMAD paper.
 {SIGMAs.push_back(0.5*tilde_Ss[0]);
  SIGMAs.push_back(0.5*tilde_Ss[1]);
  SIGMAs.push_back(0.5*tilde_Ss[2]);
  SIGMAs.push_back(0.5*tilde_Ss[3]);
  SIGMAs.push_back(DSs[0]);
  SIGMAs.push_back(DSs[1]);
  SIGMAs.push_back(DSs[2]);
  SIGMAs.push_back(DSs[3]);
  SIGMAs.push_back(DSs[4]);
  SIGMAs.push_back(DSs[5]);}

  m4d J;set_J(J); // xi matrix of the SMAD paper.
  std::vector<m4d> tilde_Ss_J(tilde_Ss); //4
  inlib::multiply(tilde_Ss_J,J);
  std::vector<m4d> DSs_J(DSs); //6
  inlib::multiply(DSs_J,J);
  std::vector<m4d> SIGMAs_J(SIGMAs); //10
  inlib::multiply(SIGMAs_J,J);
  
  //////////////////////////////////////////////////////////
  /// D=5 material: ////////////////////////////////////////
  //////////////////////////////////////////////////////////
  typedef inlib::mat<double,5> m5d;
  std::vector<m5d> Es_5;set_Es_5(Es_5);  //10 5x5 matrices.

  m5d eta_131;inlib::set_eta_131(eta_131);
  std::vector<m5d> Es_5_eta_131;  //10 5x5
 {for(unsigned int A=0;A<10;A++) Es_5_eta_131.push_back(Es_5[A]*eta_131);}

  //////////////////////////////////////////////////////////
  /// D=6 material: ////////////////////////////////////////
  //////////////////////////////////////////////////////////
  typedef inlib::mat<double,6> m6d;
  m6d I6;I6.set_identity();
  double lor_csts[6][6][6]; //The Lorentz group constants.
  inlib::set_Es_eta_group_constants(lor_csts);
  std::vector<m6d> lor_adjoint;
  inlib::set_from_csts_down<6, std::vector<m6d> >(lor_adjoint,lor_csts);
  m6d H6;inlib::set_H6(H6);
  m6d G6;set_G6(G6);

  //////////////////////////////////////////////////////////
  /// D=10 material: ///////////////////////////////////////
  //////////////////////////////////////////////////////////
  typedef inlib::mat<double,10> m10d;
  m10d eta_G6;inlib::set_eta_G6(eta_G6);
  
  //////////////////////////////////////////////////////////
  /// checks on the Es*eta matrices : //////////////////////
  //////////////////////////////////////////////////////////
  if(!inlib::check_rep_commutator<m4d,6>(std::cout,Es_eta,lor_csts,inlib::dfabs)) return EXIT_FAILURE;
  m6d minus_two_I6(I6);
  minus_two_I6.multiply(-2.0);
  if(!inlib::check_rep_metric<m4d,m6d,double>(std::cout,Es,minus_two_I6,inlib::dfabs)) return EXIT_FAILURE;

  if(!check_exp_LRC_relation_down(std::cout,Es_eta,Es_eta,minus_Es_eta,lor_adjoint,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;
  
  typedef inlib::mat<z,4> m4z;
  //////////////////////////////////////////////////////////
  /// checks on Dirac Majorana down matrices: //////////////
  //////////////////////////////////////////////////////////
  std::vector<m4z> down_maj_diracs; //4 matrices.
  inlib::set_majoranas_down(down_maj_diracs);
  if(verbose) inlib::mats_dump(std::cout,"Dirac Majorana down matrices:",down_maj_diracs);
  if(!check_anticometric(std::cout,down_maj_diracs,z(1,0),eta,inlib::cfabs)) return EXIT_FAILURE;
  
  std::vector<m4z> Ds; //6 matrices. Lorentz-spinor representation.
  inlib::set_rot_boost_A_rep(down_maj_diracs,Es,1.0/8.0,Ds);
//if(verbose) inlib::mats_dump(std::cout,"Ds with dirac Majorana down matrices:",Ds);
  
  if(!inlib::check_rep_commutator<m4z,6>(std::cout,Ds,lor_csts,inlib::cfabs)) return EXIT_FAILURE;
  if(!inlib::check_rep_gam_coms_down<m4z,m4d>(std::cout,Ds,down_maj_diracs,Es_eta,inlib::cfabs)) return EXIT_FAILURE;
  std::vector<m4z> minus_Ds(Ds);
  inlib::multiply(minus_Ds,z(-1));
  if(!check_exp_LRC_relation_down(std::cout,Ds,down_maj_diracs,minus_Ds,Es_eta,rg,inlib::cfabs,40,depsilon)) return EXIT_FAILURE;
  
  if(!check_exp_LRC_relation_down(std::cout,Ds,Ds,minus_Ds,lor_adjoint,rg,inlib::cfabs,40,depsilon)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  /// checks on Dirac Majorana up matrices: ////////////////
  //////////////////////////////////////////////////////////
  std::vector<m4z> up_maj_diracs; //4 matrices. Same as Itzykson-Zuber p694.
  inlib::set_majoranas_up(up_maj_diracs);
  if(verbose) inlib::mats_dump(std::cout,"Dirac Majorana up matrices:",up_maj_diracs);
  if(!check_anticometric(std::cout,up_maj_diracs,z(1,0),eta,inlib::cfabs)) return EXIT_FAILURE;
 {std::vector<m4d> eta_Es_eta;
  set_Es(eta_Es_eta);
 {for(size_t a=0;a<6;a++) eta_Es_eta[a] = eta*Es[a]*eta;}
  std::vector<m4z> Ds_up;
  inlib::set_rot_boost_A_rep(up_maj_diracs,eta_Es_eta,1.0/8.0,Ds_up);
  if(!check_equals_prec(std::cout,Ds_up,Ds,inlib::cfabs)) return EXIT_FAILURE;
  if(!inlib::check_rep_gam_coms_up<m4z,m4d>(std::cout,Ds,up_maj_diracs,minus_Es_eta,inlib::cfabs)) return EXIT_FAILURE;
  std::vector<m4z> minus_Ds_up(Ds_up);
  inlib::multiply(minus_Ds_up,z(-1));
  if(!check_exp_LRC_relation_up(std::cout,Ds_up,up_maj_diracs,minus_Ds_up,minus_Es_eta,rg,inlib::cfabs,40,depsilon)) return EXIT_FAILURE;}

  //////////////////////////////////////////////////////////
  /// checks on the tilds_Ss matrices : ////////////////////
  //////////////////////////////////////////////////////////
  if(verbose) mats_dump(std::cout,"tilde_Ss:",tilde_Ss);
  if(!check_is_symmetric(std::cout,tilde_Ss)) return EXIT_FAILURE;
 {std::vector<m4z> tmp(4);
 {for(size_t mu=0;mu<4;mu++) {inlib::to_complex(tilde_Ss_J[mu],tmp[mu]);tmp[mu] *= -zi;}}
  if(!check_equals_prec(std::cout,down_maj_diracs,tmp,inlib::cfabs)) return EXIT_FAILURE;}
  if(!check_anticometric(std::cout,tilde_Ss_J,-1.0,eta,inlib::dfabs)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  /// checks on the Ss matrices : //////////////////////////
  //////////////////////////////////////////////////////////
  if(verbose) mats_dump(std::cout,"Ss:",DSs);
  if(!check_is_symmetric(std::cout,DSs)) return EXIT_FAILURE;
 {std::vector<m4z> tmp(6);
 {for(size_t a=0;a<6;a++) {inlib::to_complex(DSs_J[a],tmp[a]);}}
  if(!check_equals_prec(std::cout,Ds,tmp,inlib::cfabs)) return EXIT_FAILURE;}
  if(!inlib::check_rep_commutator<m4d,6>(std::cout,DSs_J,lor_csts,inlib::dfabs)) return EXIT_FAILURE;

  std::vector<m4d> minus_DSs_J(DSs_J);
  inlib::multiply(minus_DSs_J,-1.0);
  if(!check_exp_LRC_relation_down(std::cout,DSs_J,tilde_Ss_J,minus_DSs_J,Es_eta,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  /// checks on the SIGMAs matrices : //////////////////////
  //////////////////////////////////////////////////////////
  if(!check_is_symmetric(std::cout,SIGMAs)) return EXIT_FAILURE;
  double SIGMAs_J_csts[10][10][10];
  if(!inlib::set_commutator_constants_prec< std::vector<m4d>, m10d, 10, double>(SIGMAs_J,SIGMAs_J_csts,depsilon,inlib::dfabs)) {
    std::cout << "SIGMAs_J_csts: set_commutator_constants failed." << std::endl;
    return EXIT_FAILURE;
  }
  if(!inlib::check_rep_commutator<m4d,10>(std::cout,SIGMAs_J,SIGMAs_J_csts,inlib::dfabs)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  /// checks D=5 : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
  double Es_5_eta_131_csts[10][10][10];
  if(!inlib::set_commutator_constants_prec< std::vector<m5d>, m10d, 10, double>(Es_5_eta_131,Es_5_eta_131_csts,depsilon,inlib::dfabs)) {
    std::cout << "Es_5_eta_131_csts: set_commutator_constants failed." << std::endl;
    return EXIT_FAILURE;
  }
  if(!inlib::check_rep_commutator<m5d,10>(std::cout,Es_5_eta_131,Es_5_eta_131_csts,inlib::dfabs)) return EXIT_FAILURE;

  if(!inlib::check_constants_are_equal<10,double>(std::cout,Es_5_eta_131_csts,SIGMAs_J_csts)) return EXIT_FAILURE;

  std::vector<m5d> Es_5_eta_131_4;  //4 5x5
 {for(unsigned int mu=0;mu<4;mu++) Es_5_eta_131_4.push_back(Es_5_eta_131[mu]);}
  std::vector<m5d> Es_5_eta_131_6;  //6 5x5
 {for(unsigned int a=0;a<6;a++) Es_5_eta_131_6.push_back(Es_5_eta_131[4+a]);}
  if(!inlib::check_rep_commutator<m5d,6>(std::cout,Es_5_eta_131_6,lor_csts,inlib::dfabs)) return EXIT_FAILURE;
  if(!inlib::check_rep_gam_coms_down<m5d,m4d>(std::cout,Es_5_eta_131_6,Es_5_eta_131_4,Es_eta,inlib::dfabs)) return EXIT_FAILURE;

  std::vector<m5d> minus_Es_5_eta_131_6(Es_5_eta_131_6);
  inlib::multiply(minus_Es_5_eta_131_6,-1.0);
  if(!check_exp_LRC_relation_down(std::cout,Es_5_eta_131_6,Es_5_eta_131_4,minus_Es_5_eta_131_6,Es_eta,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  /// checks D=6 : /////////////////////////////////////////
  //////////////////////////////////////////////////////////
 {for(unsigned int j=0;j<3;j++) {
   {for(unsigned int k=0;k<3;k++) {
      for(unsigned int l=0;l<3;l++) {
        if(lor_adjoint[j].value(k,l)!=-rs[j].value(k,l)) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j].value(k+3,l+3)!=-rs[j].value(k,l)) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j].value(k,l+3)!=0) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j].value(k+3,l)!=0) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }

        if(lor_adjoint[j+3].value(k,l+3)!=rs[j].value(k,l)) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j+3].value(k+3,l)!=-rs[j].value(k,l)) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j+3].value(k,l)!=0) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
        if(lor_adjoint[j+3].value(k+3,l+3)!=0) {
          std::cout << "check lor_adjoint failed." << std::endl;
          return EXIT_FAILURE;
        }
      }
    }}
  }}

  if(!inlib::check_exp_rep_and_metric(std::cout,lor_adjoint,G6,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;
  if(!inlib::check_exp_rep_and_metric(std::cout,lor_adjoint,H6,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;
  
  std::vector<m10d> SIGMAs_J_Cs;
  inlib::set_from_csts_down<10, std::vector<m10d> >(SIGMAs_J_Cs,SIGMAs_J_csts);
  if(!inlib::check_rep_commutator<m10d,10>(std::cout,SIGMAs_J_Cs,SIGMAs_J_csts,inlib::dfabs)) return EXIT_FAILURE;
  if(!inlib::check_exp_rep_and_metric(std::cout,SIGMAs_J_Cs,eta_G6,rg,inlib::dfabs,40,depsilon)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  std::cout << "exit ..." << std::endl;

#ifdef INLIB_MEM
  }inlib::mem::balance(std::cout);
#endif //INLIB_MEM
  
  return EXIT_SUCCESS;
}
