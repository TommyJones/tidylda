
# pragma once

#include <vector>
#include <cmath>

// makes a list of indices to loop over in parallel
// use this to divide documents into batches
std::vector<std::vector<std::size_t>> allocate_batch_indices(
    const std::size_t& threads,
    const std::size_t& Nd
) {
  
  // initialize the output
  std::vector<std::vector<std::size_t>> batch_indices(threads);
  
  // initialize step size. 
  // This is a ceiling to make sure we don't accidentally drop any values
  // between 0 and Nd - 1
  std::size_t batch_size = std::ceil(
    static_cast<double>(Nd) / static_cast<double>(threads)
  ); 
  
  // start index of each batch
  std::size_t batches_start = 0;
  
  
  for (auto j = 0; j < threads; j++) {
    // stop index of each batch
    auto batches_stop = std::min(batches_start + batch_size - 1, Nd - 1);
    
    // temporary vector to hold results
    std::vector<std::size_t> tmp;
    
    // add in the correct indices
    for (std::size_t k = batches_start; k < batches_stop + 1; k++) {
      tmp.push_back (k);
    }
    
    // update the start index for the next loop
    batches_start += batch_size;
    
    // put the temporary vector in our output container
    batch_indices[j] = tmp;
    
  }
  
  return batch_indices;
}

// add a list of long integer vectors together
// use this to help collect batch-local Ck's into a global Ck
std::vector<long> add_integer_vectors(
    const std::vector<std::vector<long>>& vector_list
) {
  
  // note for performance, I am *not* checking that all vectors in vector_list
  // are of the same type and the same length. Caveat emptor!
  
  auto vec_length = vector_list[0].size();
  
  // initialize result with all values equal to zero
  std::vector<long> out(vec_length);
  
  for (auto j = 0; j < vec_length; j++) {
    out[j] = 0;
  }
  
  for (auto i = 0; i < vector_list.size(); i++) {
    for (auto j = 0; j < vec_length; j++) {
      out[j] += vector_list[i][j];
    }
  }
  
  return out;
  
}

// add a list of long integer "matrices" together
// in this case, a "matrix" is a vector of vectors
std::vector<std::vector<long>> add_integer_matrices( // 
    const std::vector<std::vector<std::vector<long>>>& matrix_list
) {
  
  // note for performance, I am *not* checking format of inputs. Caveat emptor!
  
  // initialize some variables
  auto num_mats = matrix_list.size();
  
  // get "matrix" dimensions
  auto ncol = matrix_list[0].size();
  
  // make a temporary container where we will put the jth column of each matrix
  std::vector<std::vector<long>> tmp1(num_mats);
  
  // initialize output
  std::vector<std::vector<long>> out(ncol);
  
  // loop to fill her out
  for (auto k = 0; k < ncol; k++) {
    // populate tmp1 with the kth column of each matrix
    for (auto j = 0; j < num_mats; j++) {
      tmp1[j] = matrix_list[j][k];
    }
    // add them together
    out[k] = add_integer_vectors(tmp1);
  }
  
  // return the result
  return out;
}


// update global Ck from batch Ck's
std::vector<long> update_global_Ck(
  const std::vector<long>&              Ck,
  const std::vector<std::vector<long>>& Ck_batch,
  const std::size_t& threads
) {
  
  // initialize output vector 
  std::vector<long> out(Ck.size());
  
  // add together all vectors in Ck_batch
  std::vector<long> Ck_batch_sum = add_integer_vectors(Ck_batch);
  
  // calculate the output
  for (auto j = 0; j < Ck.size(); j++) {
    out[j] = Ck_batch_sum[j] - 
      static_cast<long>(threads) * Ck[j] +
      Ck[j];
  }
  
  return out;
}

// update global Cv from batch Cv's
std::vector<std::vector<long>> update_global_Cv(
    const std::vector<std::vector<long>>&              Cv,
    const std::vector<std::vector<std::vector<long>>>& Cv_batch,
    const std::size_t& threads
) {
  
  // initialize output vector 
  std::vector<std::vector<long>> out(Cv.size());
  
  // add together all "matrices" in Cv_batch
  std::vector<std::vector<long>> Cv_batch_sum = add_integer_matrices(Cv_batch);
  
  // calculate the output
  for (auto j = 0; j < Cv.size(); j++) {
    
    std::vector<long> tmp(Cv[j].size());
    
    for (auto k = 0; k < tmp.size(); k++) {
      tmp[k] = Cv_batch_sum[j][k] - 
        static_cast<long>(threads) * Cv[j][k] +
        Cv[j][k];
    }
    
    out[j] = tmp;
  }
  
  return out;
}
