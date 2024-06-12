#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mpi.h>
#include <gtest/gtest.h>

#include "../include/privateer/privateer.hpp"
#include "../test_apps/utility/random.hpp"
/*
  TODO:
  - fix incremental snapshot test
  - check r/w only permissions, create/open
  - test simple utility functions, return values of create/open
  - parameterize necessary values (the ones with 10)
  - stash based constructor
*/
std::vector<size_t> get_random_offsets(size_t region_length, size_t num_updates){
  std::vector<size_t> random_values;
  std::generate_n(std::back_inserter(random_values), num_updates, utility::RandomNumberBetween(0,region_length - 1));
  return random_values;
}

std::vector<size_t> get_random_offsets(size_t region_start, size_t region_end, size_t num_updates){
  std::vector<size_t> random_values;
  std::generate_n(std::back_inserter(random_values), num_updates, utility::RandomNumberBetween(region_start,region_end - 1));
  return random_values;
}

class PrivateerTest : public testing::TestWithParam<std::tuple<size_t, size_t, size_t, size_t>> {
  public:
    Privateer* priv = nullptr;
    size_t size_bytes;
    size_t num_ints;
    size_t* data;

    void SetUp() override {
      char env[] = "PRIVATEER_MAX_MEM_BLOCKS=";
      char block_num[10];
      strcpy(block_num, (std::to_string(std::get<0>(GetParam()))).c_str());
      strcat(env, block_num);
      putenv(env);

      priv = new Privateer(Privateer::CREATE, "datastore");
      size_bytes = std::get<1>(GetParam());
      num_ints = size_bytes / sizeof(size_t);
      data = (size_t*) priv->create(nullptr, "v0", size_bytes);
    }

    void TearDown() override {
      delete priv;
      std::filesystem::remove_all("datastore");
    }
};

// TO BE REMOVED
TEST_P(PrivateerTest, Immutable) {
  size_t start = 0;
  size_t middle = this->num_ints / 2;
  size_t middle_to_end = ( this->num_ints / 2 ) + ( this->num_ints / 4 );
  size_t end = this->num_ints - 1;
  this->data[start] = 7;
  this->data[middle] = 8;
  this->data[middle_to_end] = 9;
  this->data[end] = 10;
  std::cout << "written to: " << end << std::endl;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_immutable(nullptr, "v0", "v1");
  this->data[start] = 1;
  this->data[middle] = 2;
  this->data[middle_to_end] = 3;
  this->data[end] = 4;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");
  EXPECT_EQ(this->data[start], 7);
  EXPECT_EQ(this->data[middle], 8);
  EXPECT_EQ(this->data[middle_to_end], 9);
  EXPECT_EQ(this->data[end], 10);
  this->data = (size_t*) priv->open_read_only(nullptr, "v1");
  EXPECT_EQ(this->data[start], 1);
  EXPECT_EQ(this->data[middle], 2);
  EXPECT_EQ(this->data[middle_to_end], 3);
  EXPECT_EQ(this->data[end], 4);
}

TEST_P(PrivateerTest, SimpleWrite) {
  size_t start = 0;
  size_t middle = this->num_ints / 2;
  size_t middle_to_end = ( this->num_ints / 2 ) + ( this->num_ints / 4 );
  size_t end = this->num_ints - 1;
  this->data[start] = 7;
  this->data[middle] = 8;
  this->data[middle_to_end] = 9;
  this->data[end] = 10;
  std::cout << "written to: " << end << std::endl;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");
  EXPECT_EQ(this->data[start], 7);
  EXPECT_EQ(this->data[middle], 8);
  EXPECT_EQ(this->data[middle_to_end], 9);
  EXPECT_EQ(this->data[end], 10);
  std::cout << "region size: " << this->priv->region_size() << std::endl;
  std::cout << "version capacity: " << this->priv->version_capacity("datastore/v0") << std::endl;
  std::cout << "version block size: " << this->priv->version_block_size("datastore/v0") << std::endl;
}

TEST_P(PrivateerTest, SimpleDenseWrite) {
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = i;
  }
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");
  for (size_t i = 0; i < this->num_ints; i++) {
    EXPECT_EQ(this->data[i], i);
  }
}

TEST_P(PrivateerTest, SimpleWrite_Data) {
  size_t start = 0;
  size_t middle = this->num_ints / 2;
  size_t middle_to_end = ( this->num_ints / 2 ) + ( this->num_ints / 4 );
  size_t end = this->num_ints - 1;
  this->data[start] = 7;
  this->data[middle] = 8;
  EXPECT_EQ(this->data[start], 7);
  EXPECT_EQ(this->data[middle], 8);
  EXPECT_EQ(this->data[middle_to_end], 0);
  EXPECT_EQ(this->data[end], 0);
  EXPECT_EQ(((size_t*) this->priv->data())[start], 7);
  EXPECT_EQ(((size_t*) this->priv->data())[middle], 8);
  EXPECT_EQ(((size_t*) this->priv->data())[middle_to_end], 0);
  EXPECT_EQ(((size_t*) this->priv->data())[end], 0);
  this->data[middle_to_end] = 9;
  this->data[end] = 10;
  std::cout << "written to: " << end << std::endl;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");
  EXPECT_EQ(this->data[start], 7);
  EXPECT_EQ(this->data[middle], 8);
  EXPECT_EQ(this->data[middle_to_end], 9);
  EXPECT_EQ(this->data[end], 10);
  EXPECT_EQ(((size_t*) this->priv->data())[start], 7);
  EXPECT_EQ(((size_t*) this->priv->data())[middle], 8);
  EXPECT_EQ(((size_t*) this->priv->data())[middle_to_end], 9);
  EXPECT_EQ(((size_t*) this->priv->data())[end], 10);
}

TEST_P(PrivateerTest, SortWrite) {
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = (this->num_ints - 1) - i;
  }
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open  (nullptr, "v0");
  for (size_t i = 0; i < this->num_ints; i++) {
    EXPECT_EQ(this->data[i], (this->num_ints - 1) - i);
  }
  std::sort(this->data, this->data + this->num_ints);
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");
  for (size_t i = 0; i < this->num_ints; i++) {
    EXPECT_EQ(this->data[i], i);
  }
}

TEST_P(PrivateerTest, MultipleWrite) {
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = i;
  }

  int num_iterations = std::get<2>(GetParam());
  for (size_t k = 0; k < num_iterations; k++) {
    delete priv;
    priv = new Privateer(Privateer::OPEN, "datastore");
    this->data = (size_t*) priv->open_read_only(nullptr, "v0");
    for (size_t i = 0; i < this->num_ints; i++) {
      EXPECT_EQ(this->data[i], i);
    }
    for (size_t i = 0; i < this->num_ints; i++) {
      this->data[i] = (this->num_ints - 1) - i;
    }
    priv->msync();
    delete priv;

    priv = new Privateer(Privateer::OPEN, "datastore");
    this->data = (size_t*) priv->open  (nullptr, "v0");
    for (size_t i = 0; i < this->num_ints; i++) {
      EXPECT_EQ(this->data[i], (this->num_ints - 1) - i);
    }
    for (size_t i = 0; i < this->num_ints; i++) {
      this->data[i] = i;
    }
    priv->msync();
  }
}

TEST_P(PrivateerTest, IncrementalRandomSparseWrite) {
  int num_iterations = std::get<2>(GetParam());
  int num_updates = std::get<3>(GetParam());

  for (int i = 0 ; i < num_iterations ; i++) {
    std::vector<size_t> offsets = get_random_offsets(this->num_ints, num_updates);
    for (auto offset : offsets) {
    EXPECT_GE(offset, 0);
    EXPECT_LT(offset, this->num_ints);
      data[offset] += 1;
    }
    priv->msync();
  }
}

TEST_P(PrivateerTest, IncrementalRandomSparseWrite_Threaded) {
  int num_iterations = std::get<2>(GetParam());
  int num_updates = std::get<3>(GetParam());
  omp_set_num_threads(4);
  for (int i = 0 ; i < num_iterations ; i++) {
    std::vector<size_t> offsets = get_random_offsets(this->num_ints, num_updates);
    std::vector<size_t>::iterator offset_iterator;
#pragma omp for
    for (offset_iterator = offsets.begin(); offset_iterator < offsets.end(); ++offset_iterator){
      std::cout << "offset iterator: " << *offset_iterator << std::endl;
      EXPECT_GE(*offset_iterator, 0);
      EXPECT_LT(*offset_iterator, this->num_ints);
      data[*offset_iterator] += 1;
    }
    priv->msync();
  }
}

TEST_P(PrivateerTest, SimpleSnapshot) {
  int num_iterations = std::get<2>(GetParam());
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = 0;
  }
  priv->msync();
  for (int j = 1; j <= num_iterations; ++j){
    for (size_t k = 1; k < this->num_ints; k+=2) {
      this->data[k]++;
    }
    priv->snapshot(("v" + std::to_string(j)).c_str());
  }
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  for (int j = 1; j <= num_iterations; ++j){
    this->data = (size_t*) priv->open_read_only(nullptr, ("v" + std::to_string(j)).c_str());
    for (size_t k = 1; k < this->num_ints; k+=2){
      EXPECT_EQ(data[k], j);
      EXPECT_EQ(data[k-1], 0);
    }
  }
}

TEST_P(PrivateerTest, IncrementalRandomSparseSnapshot) {
  int num_iterations = std::get<2>(GetParam());
  size_t update_ratio = std::get<3>(GetParam());

  float initial_fill_ratio = 0.01;
  size_t initial_fill_size = this->num_ints * initial_fill_ratio;
  float initial_sparsity = 0.01;
  size_t num_updates = initial_fill_size*initial_sparsity;

  std::vector<size_t> random_indices_first_half = get_random_offsets(0, initial_fill_size, num_updates);
  for (auto offset_iterator : random_indices_first_half) {
    EXPECT_GE(offset_iterator, 0);
    EXPECT_LT(offset_iterator, this->num_ints);
    this->data[offset_iterator] += 1;
  }

  size_t update_size = num_ints*(update_ratio * 1.0 / 100);

  for (int i = 1; i < num_iterations; i++) {
  std::cout << "iteration: " << i << std::endl;
    size_t update_start = initial_fill_size + i*update_size;
    num_updates = update_size*initial_sparsity;

    std::vector<size_t> random_indices = get_random_offsets(update_start, update_start + update_size, num_updates);
    for (auto offset_iterator : random_indices) {
    EXPECT_GE(offset_iterator, 0);
    EXPECT_LT(offset_iterator, this->num_ints);
      this->data[offset_iterator] += 1;
    }

    EXPECT_TRUE(priv->snapshot(("v" + std::to_string(i)).c_str()));
  }
}

TEST_P(PrivateerTest, IncrementalRandomSparseSnapshot_Threaded) {
  int num_iterations = std::get<2>(GetParam());
  size_t update_ratio = std::get<3>(GetParam());

  float initial_fill_ratio = 0.01;
  size_t initial_fill_size = this->num_ints * initial_fill_ratio;
  float initial_sparsity = 0.01;
  size_t num_updates = initial_fill_size*initial_sparsity;

  std::vector<size_t> random_indices_first_half = get_random_offsets(0, initial_fill_size, num_updates);
  std::vector<size_t>::iterator offset_iterator;
  #pragma omp parallel for
  for (offset_iterator = random_indices_first_half.begin(); offset_iterator <= random_indices_first_half.end(); ++offset_iterator){
    EXPECT_GE(*offset_iterator, 0);
    EXPECT_LT(*offset_iterator, this->num_ints);
    std::cout << "offset_iterator: " << *offset_iterator << std::endl;
    this->data[*offset_iterator] += 1;
  }

  size_t update_size = num_ints*(update_ratio * 1.0 / 100);

  for (int i = 1; i < num_iterations; i++) {
    size_t update_start = initial_fill_size + i*update_size;
    num_updates = update_size*initial_sparsity;

    std::vector<size_t> random_indices = get_random_offsets(update_start, update_start + update_size, num_updates);
    #pragma omp parallel for
    for (offset_iterator = random_indices.begin(); offset_iterator < random_indices.end(); ++offset_iterator){
      EXPECT_GE(*offset_iterator, 0);
      EXPECT_LT(*offset_iterator, this->num_ints);
      std::cout << "offset_iterator: " << *offset_iterator << std::endl;
      this->data[*offset_iterator] += 1;
    }

    EXPECT_TRUE(priv->snapshot(("v" + std::to_string(i)).c_str()));
  }
}

TEST(PrivateerTest_Concurrent, ConcurrentWrite) {
  size_t size_bytes = 1024LLU;
  size_t num_ints = size_bytes / sizeof(size_t);
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0){
      Privateer privateer(Privateer::CREATE, "datastore");
      privateer.create(nullptr, "v0", size_bytes);
      privateer.msync();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank != 0){
      Privateer privateer(Privateer::OPEN, "datastore");
      size_t *data = (size_t*) privateer.open_immutable(nullptr, "v0", ("v" + std::to_string(world_rank)).c_str());
      for (size_t i = 0; i < num_ints; i++){
        data[i] = i;
      }
      privateer.msync();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0) {
    for (int i = 1; i < world_size; i++){
      Privateer privateer(Privateer::OPEN, "datastore");
      size_t *data = (size_t*) privateer.open_read_only(nullptr, ("v" + std::to_string(i)).c_str());
      for (size_t j = 0; j < num_ints; j++){
        EXPECT_EQ(data[j], j);
      }
    }
  }
  MPI_Finalize();
  std::filesystem::remove_all("datastore");
}

// Death tests
TEST_P(PrivateerTest, LowerBoundOutOfRange) {
  EXPECT_DEATH({
      this->data[-1] = 1;
    }, "Fault address out of range");
}

TEST_P(PrivateerTest, UpperBoundOutOfRange) {
  std::cout << "num_ints: " << this->num_ints << std::endl;
  std::cout << "size_bytes: " << this->size_bytes << std::endl;
  std::cout << "size_t size: " << sizeof(size_t) << std::endl;
  EXPECT_DEATH({
      this->data[this->num_ints] = 1;
    }, "Fault address out of range");
}

TEST_P(PrivateerTest, ReadOnly) {
  this->data[0] = 7;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  this->data = (size_t*) priv->open_read_only(nullptr, "v0");

  EXPECT_DEATH({
      this->data[0] = 1;
    }, "");
}

#ifdef USE_COMPRESSION
TEST_P(PrivateerTest, SimpleCompressionTest) {
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = 7;
  }
  priv->msync();
  size_t size = 0;
  for (const auto& entry : std::filesystem::recursive_directory_iterator("datastore/blocks")) {
    if (std::filesystem::is_regular_file(entry.path())) {
      size += std::filesystem::file_size(entry.path());
    }
  }
  EXPECT_LT(size, 2097152);
}
#endif

/*
** PARAMS
** 0 - max number of 2MB blocks, constraining max allotted memory region
** 1 - size of datastore region
** 2 - iterations
** 3 - update rate/ratio
*/

INSTANTIATE_TEST_SUITE_P(
    Parameterized_PrivateerTest,
    PrivateerTest,
    ::testing::Values(
      std::make_tuple(    1,               8 * 1024LLU, 10, 10),
      std::make_tuple(    2,               8 * 1024LLU, 10, 10),
      std::make_tuple(    4,               8 * 1024LLU, 10, 10),
      std::make_tuple(    8,               8 * 1024LLU, 10, 10),
      std::make_tuple(   16,               8 * 1024LLU, 10, 10),
      std::make_tuple(16384,               8 * 1024LLU, 10, 10),
      std::make_tuple(    1,        8 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    2,        8 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    4,        8 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    8,        8 * 1024 * 1024LLU, 10, 10),
      std::make_tuple(   16,        8 * 1024 * 1024LLU, 10, 10),
      std::make_tuple(16384,        8 * 1024 * 1024LLU, 10, 10)/*,
      std::make_tuple(    1, 8 * 1024 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    2, 8 * 1024 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    4, 8 * 1024 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(    8, 8 * 1024 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(   16, 8 * 1024 * 1024 * 1024LLU, 10, 10), // page eviction occurs
      std::make_tuple(16384, 8 * 1024 * 1024 * 1024LLU, 10, 10)*/
    )
);
