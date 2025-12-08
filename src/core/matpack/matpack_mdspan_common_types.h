#pragma once

#include <configtypes.h>

#include <experimental/mdspan>
#include <ranges>
#include <type_traits>
#include <utility>

#include "matpack_mdspan_common_sizes.h"

namespace stdx = std::experimental;
namespace stdr = std::ranges;
namespace stdv = std::ranges::views;

namespace matpack {
//! A standard type holding a strided multidimensional array
template <class T, Size N>
using mdstrided_t =
    stdx::mdspan<T, stdx::dextents<Index, N>, stdx::layout_stride>;

//! A standard type holding a multidimensional array
template <class T, Size N>
using mdview_t = stdx::mdspan<T, stdx::dextents<Index, N>>;

//! Our data holder
template <typename T, Size N>
class data_t;

//! Our view holder
template <typename T, Size N>
struct view_t;

//! Our strided view holder
template <typename T, Size N>
struct strided_view_t;

//! Our constant-sized data holder
template <typename T, Size... N>
struct cdata_t;

//! Helper traits
namespace {
template <typename T, Size... dims>
consteval std::true_type is_data_t_helper(const data_t<T, dims...>*);
consteval std::false_type is_data_t_helper(...);
template <typename T, Size... dims>
consteval std::true_type is_view_t_helper(const view_t<T, dims...>*);
consteval std::false_type is_view_t_helper(...);
template <typename T, Size... dims>
consteval std::true_type is_strided_view_t_helper(
    const strided_view_t<T, dims...>*);
consteval std::false_type is_strided_view_t_helper(...);
template <typename T, Size... dims>
consteval std::true_type is_cdata_t_helper(const cdata_t<T, dims...>*);
consteval std::false_type is_cdata_t_helper(...);
}  // namespace

//! Type trait to detect cdata_t or types derived from data_t
template <typename T>
struct is_data_t
    : decltype(is_data_t_helper(std::declval<std::remove_cvref_t<T>*>())){};

//! Type trait to detect strided_view_t or types derived from view_t
template <typename T>
struct is_view_t
    : decltype(is_view_t_helper(std::declval<std::remove_cvref_t<T>*>())){};

//! Type trait to detect strided_view_t or types derived from strided_view_t
template <typename T>
struct is_strided_view_t : decltype(is_strided_view_t_helper(
                               std::declval<std::remove_cvref_t<T>*>())){};

//! Type trait to detect cdata_t or types derived from cdata_t
template <typename T>
struct is_cdata_t
    : decltype(is_cdata_t_helper(std::declval<std::remove_cvref_t<T>*>())){};

template <typename T>
concept any_data = is_data_t<std::remove_cvref_t<T>>::value;

template <typename T>
concept any_view = is_view_t<std::remove_cvref_t<T>>::value;

template <typename T>
concept any_strided_view = is_strided_view_t<std::remove_cvref_t<T>>::value;

template <typename T>
concept any_cdata = is_cdata_t<std::remove_cvref_t<T>>::value;

//! Collection of all matpack core types
template <class T>
concept any_md =
    any_data<T> or any_strided_view<T> or any_view<T> or any_cdata<T>;

template <typename T>
concept has_value_type =
    requires { typename std::remove_cvref_t<T>::value_type; };

template <has_value_type T>
using value_type = typename std::remove_cvref_t<T>::value_type;

template <typename T>
concept has_is_const = requires { std::remove_cvref_t<T>::is_const; };

template <has_is_const T>
constexpr bool mdmutable = not std::remove_cvref_t<T>::is_const;

template <typename T>
concept const_forwarded =
    std::is_const_v<T> or std::is_const_v<std::remove_reference_t<T>>;

template <typename T>
concept has_element_type =
    requires { typename std::remove_cvref_t<T>::element_type; };

template <has_element_type T>
using element_type = typename std::remove_cvref_t<T>::element_type;

template <typename T>
concept has_mapping_type =
    requires { typename std::remove_cvref_t<T>::mapping_type; };

template <has_element_type T>
using mapping_type = typename std::remove_cvref_t<T>::mapping_type;

////////////////////////////////////////////////////////////////////////////////
// Any matpack core type - /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T>
concept mut_any_md = any_md<T> and mdmutable<T>;

////////////////////////////////////////////////////////////////////////////////
// Any matpack core type of rank N - ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T, Size N>
concept ranked_data = any_data<T> and rank<T>() == N;

template <class T, Size N>
concept ranked_view = any_view<T> and rank<T>() == N;

template <class T, Size N>
concept ranked_strided_view = any_strided_view<T> and rank<T>() == N;

template <class T, Size N>
concept ranked_cdata = any_cdata<T> and rank<T>() == N;

template <typename T, Size N>
concept ranked_md = any_md<T> and rank<T>() == N;

template <typename T, Size N>
concept mut_ranked_md = ranked_md<T, N> and mdmutable<T>;

template <typename T, Size N>
concept const_ranked_md = ranked_md<T, N> and not mdmutable<T>;

template <typename T, typename U, Size N>
concept ranked_convertible_md =
    ranked_md<T, N> and std::convertible_to<value_type<T>, U>;

template <typename T, Size N>
concept dyn_ranked_md = ranked_md<T, N> and not any_cdata<T>;

////////////////////////////////////////////////////////////////////////////////
// Any matpack core type of holding unqualified U - ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T, class U>
concept typed_data = any_data<T> and std::same_as<value_type<T>, U>;

template <class T, class U>
concept typed_view = any_view<T> and std::same_as<value_type<T>, U>;

template <class T, class U>
concept typed_strided_view =
    any_strided_view<T> and std::same_as<value_type<T>, U>;

template <class T, class U>
concept typed_cdata = any_cdata<T> and std::same_as<value_type<T>, U>;

template <typename T, class U>
concept typed_md = any_md<T> and std::same_as<value_type<T>, U>;

template <typename T, typename U>
concept mut_typed_md = typed_md<T, U> and mdmutable<T>;

template <typename T, typename U>
concept const_typed_md = typed_md<T, U> and not mdmutable<T>;

////////////////////////////////////////////////////////////////////////////////
// Any matpack core type of holding unqualified U and of rank N - //////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U, Size N>
concept exact_data = ranked_data<T, N> and typed_data<T, U>;

template <typename T, typename U, Size N>
concept exact_view = ranked_view<T, N> and typed_view<T, U>;

template <typename T, typename U, Size N>
concept exact_strided_view =
    ranked_strided_view<T, N> and typed_strided_view<T, U>;

template <typename T, typename U, Size N>
concept exact_cdata = ranked_cdata<T, N> and typed_cdata<T, U>;

template <typename T, typename U, Size N>
concept exact_md = ranked_md<T, N> and typed_md<T, U>;

template <typename T, typename U, Size N>
concept mut_exact_md = exact_md<T, U, N> and mdmutable<T>;

template <typename T, typename U, Size N>
concept const_exact_md = exact_md<T, U, N> and not mdmutable<T>;
}  // namespace matpack
