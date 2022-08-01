// Glaze Library
// For the license information refer to glaze.hpp

#pragma once

#include <string>
#include <type_traits>
#include <variant>

#include "frozen/string.h"
#include "frozen/unordered_map.h"

#include "NanoRange/nanorange.hpp"

namespace glaze
{   
   namespace detail
   {
      template <int... I>
      using is = std::integer_sequence<int, I...>;
      template <int N>
      using make_is = std::make_integer_sequence<int, N>;

      constexpr auto size(const char *s) noexcept
      {
         int i = 0;
         while (*s != 0) {
            ++i;
            ++s;
         }
         return i;
      }

      template <const char *, typename, const char *, typename>
      struct concat_impl;

      template <const char *S1, int... I1, const char *S2, int... I2>
      struct concat_impl<S1, is<I1...>, S2, is<I2...>>
      {
         static constexpr const char value[]{S1[I1]..., S2[I2]..., 0};
      };

      template <const char *S1, const char *S2>
      constexpr auto concat_char()
      {
         return concat_impl<S1, make_is<size(S1)>, S2,
                            make_is<size(S2)>>::value;
      };

      template <size_t... Is>
      struct seq
      {};
      template <size_t N, size_t... Is>
      struct gen_seq : gen_seq<N - 1, N - 1, Is...>
      {};
      template <size_t... Is>
      struct gen_seq<0, Is...> : seq<Is...>
      {};

      template <size_t N1, size_t... I1, size_t N2, size_t... I2>
      constexpr std::array<char const, N1 + N2 - 1> concat(char const (&a1)[N1],
                                                           char const (&a2)[N2],
                                                           seq<I1...>,
                                                           seq<I2...>)
      {
         return {{a1[I1]..., a2[I2]...}};
      }

      template <size_t N1, size_t N2>
      constexpr std::array<char const, N1 + N2 - 1> concat_arrays(
         char const (&a1)[N1], char const (&a2)[N2])
      {
         return concat(a1, a2, gen_seq<N1 - 1>{}, gen_seq<N2>{});
      }

      template <size_t N>
      struct string_literal
      {
         static constexpr size_t size = (N > 0) ? (N - 1) : 0;
         constexpr string_literal(const char (&str)[N])
         {
            std::copy_n(str, N, value);
         }
         char value[N];
         constexpr const char *end() const noexcept { return value + size; }

         constexpr const std::string_view sv() const noexcept
         {
            return {value, size};
         }
      };

      template <size_t N>
      constexpr size_t length(char const (&)[N]) noexcept
      {
         return N;
      }

      template <string_literal Str>
      struct chars_impl
      {
         static constexpr std::string_view value{Str.value,
                                                 length(Str.value) - 1};
      };

      template <string_literal Str>
      inline constexpr std::string_view chars = chars_impl<Str>::value;
      
      template <uint32_t Format>
      struct read {};
      
      template <uint32_t Format>
      struct write {};
   }  // namespace detail

   template <class T>
   struct meta
   {};

   template <class T>
   using meta_t = std::decay_t<decltype(meta<T>::value)>;

   template <class T>
   inline constexpr auto &meta_v = meta<T>::value;

   struct comment_t
   {
      std::string_view str;
   };

   constexpr comment_t operator"" _c(const char *s, std::size_t n) noexcept
   {
      return comment_t{{s, n}};
   }

   struct raw_json
   {
      std::string str;
   };

   using basic =
      std::variant<bool, char, char8_t, unsigned char, signed char, char16_t,
                   short, unsigned short, wchar_t, char32_t, float, int,
                   unsigned int, long, unsigned long, double, long double,
                   long long, unsigned long long, std::string>;

   using basic_ptr =
      std::variant<bool *, char *, char8_t *, unsigned char *, signed char *,
                   char16_t *, short *, unsigned short *, wchar_t *, char32_t *,
                   float *, int *, unsigned int *, long *, unsigned long *,
                   double *, long double *, long long *, unsigned long long *,
                   std::string *>;

   namespace detail
   {
      template <class T>
      concept char_t = std::same_as<std::decay_t<T>, char> || std::same_as<std::decay_t<T>, char16_t> ||
         std::same_as<std::decay_t<T>, char32_t> || std::same_as<std::decay_t<T>, wchar_t>;

      template <class T>
      concept bool_t =
         std::same_as<std::decay_t<T>, bool> || std::same_as<std::decay_t<T>, std::vector<bool>::reference>;

      template <class T>
      concept int_t = std::integral<std::decay_t<T>> && !char_t<std::decay_t<T>> && !bool_t<T>;

      template <class T>
      concept num_t = std::floating_point<std::decay_t<T>> || int_t<T>;

      template <class T>
      concept glaze_t = requires
      {
         meta<T>::value;
      };

      template <class T>
      concept complex_t = glaze_t<std::decay_t<T>>;

      template <class T>
      concept str_t = !complex_t<T> && std::convertible_to<std::decay_t<T>, std::string_view>;

      template <class T>
      concept pair_t = requires(T pair)
      {
         {
            pair.first
            } -> std::same_as<typename T::first_type &>;
         {
            pair.second
            } -> std::same_as<typename T::second_type &>;
      };

      template <class T>
      concept map_subscriptable = requires(T container)
      {
         {
            container[std::declval<typename T::key_type>()]
            } -> std::same_as<typename T::mapped_type &>;
      };

      template <class T>
      concept map_t =
         !complex_t<T> && !str_t<T> && nano::ranges::range<T> &&
         pair_t<nano::ranges::range_value_t<T>> && map_subscriptable<T>;

      template <class T>
      concept array_t =
         !complex_t<T> && !str_t<T> && !map_t<T> && nano::ranges::range<T>;

      template <class T>
      concept emplace_backable = requires(T container)
      {
         {
            container.emplace_back()
            } -> std::same_as<typename T::reference>;
      };

      template <class T>
      concept resizeable = requires(T container)
      {
         container.resize(0);
      };

      template <class T>
      concept tuple_t = requires(T t)
      {
         std::tuple_size<T>::value;
         std::get<0>(t);
      }
      &&!complex_t<T> && !nano::ranges::range<T>;

      template <class T>
      concept nullable_t = !complex_t<T> && !str_t<T> && requires(T t)
      {
         bool(t);
         {*t};
      };

      template <class... T>
      constexpr bool all_member_ptr(std::tuple<T...>)
      {
         return std::conjunction_v<std::is_member_pointer<std::decay_t<T>>...>;
      }

      template <class T>
      concept glaze_array_t = glaze_t<T> && all_member_ptr(meta_v<T>);

      template <class T>
      concept glaze_object_t = glaze_t<T> && !glaze_array_t<T>;

      template <class From, class To>
      concept non_narrowing_convertable = requires(From from, To to)
      {
         To{from};
      };
      
      template <class T>
      concept stream_t = requires(T t) {
         typename T::char_type;
         typename T::traits_type;
         typename T::int_type;
         t.get();
         t.peek();
         t.unget();
         t.gcount();
      };

      // from
      // https://stackoverflow.com/questions/55941964/how-to-filter-duplicate-types-from-tuple-c
      template <class T, class... Ts>
      struct unique
      {
         using type = T;
      };

      template <template <class...> class T, class... Ts, class U, class... Us>
      struct unique<T<Ts...>, U, Us...>
         : std::conditional_t<(std::is_same_v<U, Ts> || ...),
                              unique<T<Ts...>, Us...>,
                              unique<T<Ts..., U>, Us...>>
      {};

      template <class T>
      struct tuple_variant;

      template <class... Ts>
      struct tuple_variant<std::tuple<Ts...>> : unique<std::variant<>, Ts...>
      {};

      template <class T>
      struct tuple_ptr_variant;

      template <class... Ts>
      struct tuple_ptr_variant<std::tuple<Ts...>>
         : unique<std::variant<>, std::add_pointer_t<Ts>...>
      {};

      template <class Tuple,
                class = std::make_index_sequence<std::tuple_size<Tuple>::value>>
      struct value_tuple_variant;

      template <class Tuple, size_t... I>
      struct value_tuple_variant<Tuple, std::index_sequence<I...>>
      {
         using type = typename tuple_variant<decltype(std::tuple_cat(
            std::declval<std::tuple<std::tuple_element_t<
               1, std::tuple_element_t<I, Tuple>>>>()...))>::type;
      };

      template <class Tuple>
      using value_tuple_variant_t = typename value_tuple_variant<Tuple>::type;

      // from
      // https://stackoverflow.com/questions/16337610/how-to-know-if-a-type-is-a-specialization-of-stdvector
      template <typename Test, template <typename...> class Ref>
      struct is_specialization : std::false_type
      {};

      template <template <typename...> class Ref, typename... Args>
      struct is_specialization<Ref<Args...>, Ref> : std::true_type
      {};

      template <class T, size_t... I>
      inline constexpr auto make_array_impl(std::index_sequence<I...>)
      {
         using value_t = typename tuple_variant<meta_t<T>>::type;
         return std::array<value_t, std::tuple_size_v<meta_t<T>>>{
            std::get<I>(meta_v<T>)...};
      }

      template <class T>
      inline constexpr auto make_array()
      {
         constexpr auto indices =
            std::make_index_sequence<std::tuple_size_v<meta_t<T>>>{};
         return make_array_impl<T>(indices);
      }

      template <class Tuple, std::size_t... Is>
      inline constexpr auto tuple_runtime_getter(std::index_sequence<Is...>)
      {
         using value_t = typename tuple_ptr_variant<Tuple>::type;
         using tuple_ref = std::add_lvalue_reference_t<Tuple>;
         using getter_t = value_t (*)(tuple_ref);
         return std::array<getter_t, std::tuple_size_v<Tuple>>{
            +[](tuple_ref tuple) -> value_t {
               return &std::get<Is>(tuple);
            }...};
      }

      template <class Tuple>
      inline auto get_runtime(Tuple &&tuple, const size_t index)
      {
         static constexpr auto indices =
            std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{};
         static constexpr auto runtime_getter =
            tuple_runtime_getter<std::decay_t<Tuple>>(indices);
         return runtime_getter[index](tuple);
      }

      template <class M>
      inline constexpr void check_member()
      {
         static_assert(std::tuple_size_v<M> == 0 || std::tuple_size_v<M> > 1,
                       "members need at least a name and a member pointer");
         static_assert(
            std::tuple_size_v<M> < 4,
            "only member_ptr, name, and comment are supported at the momment");
         if constexpr (std::tuple_size_v < M >> 0)
            static_assert(str_t<std::tuple_element_t<0, M>>,
                          "first element should be the name");
         if constexpr (std::tuple_size_v < M >> 1)
            static_assert(std::is_member_pointer_v<std::tuple_element_t<1, M>>,
                          "second element should be the member pointer");
         if constexpr (std::tuple_size_v < M >> 2)
            static_assert(std::is_same_v<std::tuple_element_t<2, M>, comment_t>,
                          "third element should be a comment_t");
      };

      template <size_t I = 0, class Tuple, class Members = std::tuple<>,
                class Member = std::tuple<>>
      constexpr auto group_members(Tuple &&tuple, Members &&members = {},
                                          Member &&member = {})
      {
         if constexpr (std::tuple_size_v<Tuple> == 0) {
            return std::make_tuple();
         }
         else if constexpr (I >= std::tuple_size_v<Tuple>) {
            check_member<Member>();
            return std::tuple_cat(members, std::make_tuple(member));
         }
         else if constexpr (I == 0) {
            return group_members<I + 1>(std::forward<Tuple>(tuple),
                                        std::forward<Members>(members),
                                        std::make_tuple(std::get<I>(tuple)));
         }
         else if constexpr (str_t<std::tuple_element_t<I, Tuple>>) {
            check_member<Member>();
            return group_members<I + 1>(
               std::forward<Tuple>(tuple),
               std::tuple_cat(members, std::make_tuple(member)),
               std::make_tuple(std::get<I>(tuple)));
         }
         else {
            return group_members<I + 1>(
               std::forward<Tuple>(tuple), std::forward<Members>(members),
               std::tuple_cat(member, std::make_tuple(std::get<I>(tuple))));
         }
      }

      template <class T, size_t... I>
      constexpr auto make_map_impl(std::index_sequence<I...>)
      {
         using value_t = value_tuple_variant_t<meta_t<T>>;
         return frozen::make_unordered_map<frozen::string, value_t,
                                           std::tuple_size_v<meta_t<T>>>(
            {std::make_pair<frozen::string, value_t>(
               frozen::string(std::get<0>(std::get<I>(meta_v<T>))),
               std::get<1>(std::get<I>(meta_v<T>)))...});
      }

      template <class T>
      constexpr auto make_map()
      {
         constexpr auto indices =
            std::make_index_sequence<std::tuple_size_v<meta_t<T>>>{};
         return make_map_impl<T>(indices);
      }
      
      template <class T, size_t... I>
      constexpr auto make_int_map_impl(std::index_sequence<I...>)
      {
         using value_t = value_tuple_variant_t<meta_t<T>>;
         return frozen::make_unordered_map<size_t, value_t,
                                           std::tuple_size_v<meta_t<T>>>(
            {std::make_pair<size_t, value_t>(
               I,
               std::get<1>(std::get<I>(meta_v<T>)))...});
      }
      
      template <class T>
      constexpr auto make_int_map()
      {
         constexpr auto indices =
            std::make_index_sequence<std::tuple_size_v<meta_t<T>>>{};
         return make_map_impl<T>(indices);
      }
      
      template <class T = void>
      struct from_json {};
      
      template <class T = void>
      struct to_json {};
   }  // namespace detail

   constexpr auto array(auto &&...args)
   {
      return std::make_tuple(args...);
   }

   constexpr auto object(auto&&... args)
   {
      return detail::group_members(
         std::make_tuple(args...));
   }
}  // namespace glaze