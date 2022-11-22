#pragma once

namespace iporoskun::longman::detail {

// Implementation from
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/

template<class T, class /*TagType*/>
class named_type {
public:
  explicit named_type(T const& value) noexcept
	: value_(value) {}
  explicit named_type(T&& value) noexcept
	: value_(std::move(value)) {}

  named_type() noexcept = default;
  named_type(const named_type&) noexcept = default;
  named_type(named_type&&) noexcept = default;
  named_type& operator=(named_type const&) noexcept = default;
  named_type& operator=(named_type&&) noexcept = default;
  operator T() const noexcept { return value_; };

  T& get() noexcept { return value_; }
  T const& get() const noexcept { return value_; }

private:
  T value_;
};

} // namespace iporoskun::longman::detail