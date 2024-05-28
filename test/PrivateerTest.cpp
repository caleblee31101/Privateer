#include <cstddef>
#include <filesystem>
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
class PrivateerTest : public testing::Test {
  public:
    Privateer* priv = nullptr;

    static void SetUpTestSuite() {
      std::generate_n(std::back_inserter(random_values), 100, utility::RandomNumberBetween(0,100 - 1));
    }
    static void TearDownTestSuite() {
    }
    void SetUp() override {
      priv = new Privateer(Privateer::CREATE, "datastore");
    }

    void TearDown() override {
      delete priv;
      std::filesystem::remove_all("datastore");
    }

    static std::vector<size_t> random_values;
};

std::vector<size_t> PrivateerTest::random_values;

TEST_F(PrivateerTest, SimpleWrite) {
  size_t size_bytes = 1024LLU;
  size_t num_ints = size_bytes / sizeof(size_t);
  size_t* data = (size_t*) priv->create(nullptr, "v0", size_bytes);

  size_t start = 0;
  size_t middle = num_ints / 2;
  size_t middle_to_end = ( num_ints / 2 ) + ( num_ints / 4 );
  size_t end = num_ints - 1;
  data[start] = 7;
  data[middle] = 8;
  data[middle_to_end] = 9;
  data[end] = 10;
  std::cout << "written to: " << end << std::endl;
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  data = (size_t*) priv->open_read_only(nullptr, "v0");
  EXPECT_EQ(data[start], 7);
  EXPECT_EQ(data[middle], 8);
  EXPECT_EQ(data[middle_to_end], 9);
  EXPECT_EQ(data[end], 10);
}

TEST_F(PrivateerTest, SimpleDenseWrite) {
  size_t size_bytes = 1024LLU;
  size_t num_ints = size_bytes / sizeof(size_t);
  size_t* data = (size_t*) priv->create(nullptr, "v0", size_bytes);
  for (size_t i = 0; i < num_ints; i++) {
    data[i] = i;
  }
  priv->msync();
  delete priv;

  priv = new Privateer(Privateer::OPEN, "datastore");
  data = (size_t*) priv->open_read_only(nullptr, "v0");
  for (size_t i = 0; i < num_ints; i++) {
    EXPECT_EQ(data[i], i);
  }
}

TEST_F(PrivateerTest, SparseWrite) {
  size_t size_bytes = 1024LLU;
  size_t num_ints = size_bytes / sizeof(size_t);
  size_t* data = (size_t*) priv->create(nullptr, "v0", size_bytes);
}

TEST_F(PrivateerTest, SimpleSnapshot) {

}

TEST_F(PrivateerTest, IterativeSnapshot) {

}
/*
TEST(HelloTest, basic_privateer) {
  //Privateer* priv = new Privateer(Privateer::CREATE, "datastore", "stash");
  //delete priv;

  //size_t* the_ints = (size_t*) priv.create(nullptr, "v0", 1000);
  //priv.msync();
  //std::filesystem::remove_all("datastore");
}
int main(int argc, char *argv[]) {
  std::cout << "0" << std::endl;
  Privateer* priv = new Privateer(Privateer::CREATE, "datastore");
  std::cout << "1" << std::endl;
  size_t size_bytes = 1024*1024*1024LLU;
  size_t* the_ints = (size_t*) priv->create(nullptr, "v0", size_bytes);
  std::cout << "2" << std::endl;
  priv->msync();
  std::cout << "3" << std::endl;
  delete priv;
  std::filesystem::remove_all("datastore");
  return 0;
}

*/
