#include "whipple.h"
#include "gtest/gtest.h"

class SpecificSingleton : public Singleton<SpecificSingleton>
{
  friend class Singleton<SpecificSingleton>;
 public:
  int Data() const { return Data_; }
 private:
  SpecificSingleton() : Data_(42) {}
  ~SpecificSingleton() {}
  int Data_;
};

TEST(WhippleTest, DefaultConstructor)
{
  Whipple a;
  EXPECT_EQ(a.w, 1.02);
}

