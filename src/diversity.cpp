// Copyright 2022 guillaume-gricourt
#include "diversity.hpp"

double Diversity::getMutRate() const noexcept { return mut_rate; }
double Diversity::getIndelFrac() const noexcept { return indel_frac; }
double Diversity::getIndelExtend() const noexcept { return indel_extend; }
int Diversity::getMaxInsSize() const noexcept { return max_ins_size; }

void Diversity::setMutRate(double i) { mut_rate = i; }
void Diversity::setIndelFrac(double i) { indel_frac = i; }
void Diversity::setIndelExtend(double i) { indel_extend = i; }
void Diversity::setMaxInsSize(double i) { max_ins_size = i; }
