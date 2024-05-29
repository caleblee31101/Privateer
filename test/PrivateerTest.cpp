#include <cstddef>
#include <filesystem>
#include <mpi.h>
#include <gtest/gtest.h>

#include "../include/privateer/privateer.hpp"
#include "../test_apps/utility/random.hpp"
/*
  Privateer(int action, const char *base_path);
  Privateer(int action, const char *base_path, const char *stash_base_path);
  ~Privateer();
  void *create(void *addr, const char *version_metadata_path, size_t region_size, bool allow_overwrite = true);
  void *open(void *addr, const char *version_metadata_path);
  void *open_read_only(void *addr, const char *version_metadata_path);
  void *open_immutable(void *addr, const char *version_metadata_path, const char *new_version_metadata_path);
  void msync();
  bool snapshot(const char *version_metadata_path);
  size_t get_block_size();
  void *data();
  bool version_exists(const char *version_metadata_path);
  size_t region_size();
  static size_t version_capacity(std::string version_path);
  static size_t version_block_size(std::string version_path);
*/
std::vector<size_t> get_random_offsets(size_t region_length, size_t num_updates){
  std::vector<size_t> random_values;
  std::generate_n(std::back_inserter(random_values), num_updates, utility::RandomNumberBetween(0,region_length - 1));
  return random_values;
}
class PrivateerTest : public testing::Test {
  public:
    Privateer* priv = nullptr;
    size_t size_bytes;
    size_t num_ints;
    size_t* data;

    static void SetUpTestSuite() {

    }
    static void TearDownTestSuite() {
    }
    void SetUp() override {
      priv = new Privateer(Privateer::CREATE, "datastore");
      size_bytes = 1024LLU;
      num_ints = size_bytes / sizeof(size_t);
      data = (size_t*) priv->create(nullptr, "v0", size_bytes);
    }

    void TearDown() override {
      delete priv;
      std::filesystem::remove_all("datastore");
    }
};

TEST_F(PrivateerTest, SimpleWrite) {
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
}

TEST_F(PrivateerTest, SimpleDenseWrite) {
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

TEST_F(PrivateerTest, SortWrite) {
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

TEST_F(PrivateerTest, SparseWrite) {
  int num_iterations = 10;
  int num_updates = 10;
  omp_set_num_threads(4);
  for (int i = 0 ; i < num_iterations ; i++) {
    std::vector<size_t> offsets = get_random_offsets(this->num_ints, num_updates);
    std::vector<size_t>::iterator offset_iterator;
    #pragma omp for
    for (offset_iterator = offsets.begin(); offset_iterator <= offsets.end(); ++offset_iterator){
      data[*offset_iterator] += 1;
    //for (auto offset : offsets) {
    //data[offset] += 1;
    }
    priv->msync();
  }
}

TEST_F(PrivateerTest, SimpleSnapshot) {
  for (size_t i = 0; i < this->num_ints; i++) {
    this->data[i] = 0;
  }
  priv->msync();
  for (int j = 1; j <= 10; ++j){
    for (size_t k = 1; k < this->num_ints; k+=2){
      this->data[k]++;
    }
    priv->snapshot(("v" + std::to_string(j)).c_str());
  }
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  for (int j = 1; j <= 10; ++j){
    this->data = (size_t*) priv->open_read_only(nullptr, ("v" + std::to_string(j)).c_str());
    for (size_t k = 1; k < this->num_ints; k+=2){
      EXPECT_EQ(data[k], j);
      EXPECT_EQ(data[k-1], 0);
    }
  }
}

TEST_F(PrivateerTest, IncrementalSnapshot) {

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
