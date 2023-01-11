// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_DIVERSITY_HPP_
#define INCLUDE_DIVERSITY_HPP_

class Diversity {
private:
  double mut_rate;
  double indel_frac;
  double indel_extend;
  int max_ins_size;

public:
  double getMutRate() const noexcept;
  double getIndelFrac() const noexcept;
  double getIndelExtend() const noexcept;
  int getMaxInsSize() const noexcept;

  void setMutRate(double);
  void setIndelFrac(double);
  void setIndelExtend(double);
  void setMaxInsSize(double);
};

#endif // INCLUDE_DIVERSITY_HPP_
