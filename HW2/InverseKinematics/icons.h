#pragma once
#include <cstdint>

inline constexpr int ICON_MIN = 0xf000;
inline constexpr int ICON_MAX = 0xf007;
#define ICON_PLAY "\xef\x80\x80"   // U+f000
#define ICON_PAUSE "\xef\x80\x81"  // U+f001
#define ICON_STOP "\xef\x80\x82"   // U+f002
#define ICON_PLUS "\xef\x80\x83"   // U+f003
#define ICON_MINUS "\xef\x80\x84"  // U+f004
#define ICON_EYE "\xef\x80\x85"    // U+f005
#define ICON_CODE "\xef\x80\x86"   // U+f006
// The Fork Awesome font is licensed under the SIL OFL 1.1:
//     http://scripts.sil.org/OFL
inline constexpr const uint8_t forkawesome[3004] = {
    0x00, 0x01, 0x00, 0x00, 0x00, 0x0e, 0x00, 0x80, 0x00, 0x03, 0x00, 0x60, 0x46, 0x46, 0x54, 0x4d, 0x93, 0x80, 0xd0,
    0x6c, 0x00, 0x00, 0x0b, 0xa0, 0x00, 0x00, 0x00, 0x1c, 0x47, 0x44, 0x45, 0x46, 0x00, 0x27, 0x00, 0x11, 0x00, 0x00,
    0x0b, 0x80, 0x00, 0x00, 0x00, 0x1e, 0x4f, 0x53, 0x2f, 0x32, 0x61, 0x68, 0xe9, 0xe5, 0x00, 0x00, 0x01, 0x68, 0x00,
    0x00, 0x00, 0x60, 0x63, 0x6d, 0x61, 0x70, 0xdf, 0xfa, 0x16, 0xff, 0x00, 0x00, 0x01, 0xf4, 0x00, 0x00, 0x01, 0x4a,
    0x63, 0x76, 0x74, 0x20, 0x00, 0x3b, 0x04, 0x6f, 0x00, 0x00, 0x03, 0x40, 0x00, 0x00, 0x00, 0x04, 0x67, 0x61, 0x73,
    0x70, 0xff, 0xff, 0x00, 0x03, 0x00, 0x00, 0x0b, 0x78, 0x00, 0x00, 0x00, 0x08, 0x67, 0x6c, 0x79, 0x66, 0x0f, 0x8e,
    0xbd, 0x94, 0x00, 0x00, 0x03, 0x5c, 0x00, 0x00, 0x03, 0x54, 0x68, 0x65, 0x61, 0x64, 0x1e, 0x0d, 0x59, 0xeb, 0x00,
    0x00, 0x00, 0xec, 0x00, 0x00, 0x00, 0x36, 0x68, 0x68, 0x65, 0x61, 0x0d, 0xdb, 0x05, 0x98, 0x00, 0x00, 0x01, 0x24,
    0x00, 0x00, 0x00, 0x24, 0x68, 0x6d, 0x74, 0x78, 0x30, 0x4d, 0x00, 0x3b, 0x00, 0x00, 0x01, 0xc8, 0x00, 0x00, 0x00,
    0x2c, 0x6c, 0x6f, 0x63, 0x61, 0x02, 0xe6, 0x03, 0xd0, 0x00, 0x00, 0x03, 0x44, 0x00, 0x00, 0x00, 0x18, 0x6d, 0x61,
    0x78, 0x70, 0x00, 0x50, 0x00, 0x6b, 0x00, 0x00, 0x01, 0x48, 0x00, 0x00, 0x00, 0x20, 0x6e, 0x61, 0x6d, 0x65, 0x50,
    0x0f, 0xa4, 0x38, 0x00, 0x00, 0x06, 0xb0, 0x00, 0x00, 0x04, 0x56, 0x70, 0x6f, 0x73, 0x74, 0x22, 0x23, 0x5a, 0xb3,
    0x00, 0x00, 0x0b, 0x08, 0x00, 0x00, 0x00, 0x70, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0xe0, 0x4d, 0xed,
    0xbe, 0x5f, 0x0f, 0x3c, 0xf5, 0x00, 0x1f, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0xdb, 0x18, 0xe0, 0x91, 0x00, 0x00,
    0x00, 0x00, 0xdc, 0x9b, 0x30, 0x5d, 0x00, 0x00, 0xff, 0x72, 0x07, 0x26, 0x05, 0x8e, 0x00, 0x00, 0x00, 0x08, 0x00,
    0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x06, 0x12, 0xfe, 0x66, 0x00, 0xa1, 0x07, 0x26,
    0x00, 0x00, 0x00, 0x00, 0x07, 0x26, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x3a, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00, 0x40, 0x00, 0x2e, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x05, 0x6d, 0x01, 0x90, 0x00, 0x05, 0x00, 0x00, 0x04, 0x8c, 0x04, 0xe6, 0x00, 0x00, 0x00, 0xfa, 0x04, 0x8c,
    0x04, 0xe6, 0x00, 0x00, 0x03, 0x5c, 0x00, 0x59, 0x01, 0xcf, 0x00, 0x00, 0x02, 0x00, 0x05, 0x03, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x50, 0x66, 0x45, 0x64, 0x00, 0x80, 0x00, 0x20, 0xf0, 0x06, 0x06, 0x00, 0xff, 0x00, 0x00, 0xa1, 0x06, 0x12, 0x01,
    0x9a, 0x80, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x20, 0x00, 0x01,
    0x02, 0x8b, 0x00, 0x3b, 0x00, 0x00, 0x00, 0x00, 0x02, 0x55, 0x00, 0x00, 0x00, 0xc8, 0x00, 0x00, 0x05, 0x7f, 0x00,
    0x00, 0x06, 0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x05, 0x80, 0x00, 0x00, 0x05, 0x80, 0x00, 0x00, 0x07, 0x00,
    0x00, 0x00, 0x07, 0x26, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x1c, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x44, 0x00, 0x03, 0x00, 0x01, 0x00, 0x00, 0x00, 0x1c, 0x00, 0x04, 0x00, 0x28,
    0x00, 0x00, 0x00, 0x06, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x20, 0xf0, 0x06, 0xff, 0xff, 0x00, 0x00, 0x00,
    0x20, 0xf0, 0x00, 0xff, 0xff, 0xff, 0xe3, 0x10, 0x04, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x01, 0x06, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3b, 0x04, 0x6f,
    0x00, 0x00, 0x00, 0x2c, 0x00, 0x2c, 0x00, 0x2c, 0x00, 0x2c, 0x00, 0x48, 0x00, 0x7c, 0x00, 0x9a, 0x00, 0xd0, 0x00,
    0xec, 0x01, 0x42, 0x01, 0xaa, 0x00, 0x02, 0x00, 0x3b, 0x00, 0x00, 0x02, 0x15, 0x04, 0xaa, 0x00, 0x03, 0x00, 0x07,
    0x00, 0x2e, 0xb1, 0x01, 0x00, 0x2f, 0x3c, 0xb2, 0x07, 0x04, 0x00, 0xed, 0x32, 0xb1, 0x06, 0x05, 0xdc, 0x3c, 0xb2,
    0x03, 0x02, 0x00, 0xed, 0x32, 0x00, 0xb1, 0x03, 0x00, 0x2f, 0x3c, 0xb2, 0x05, 0x04, 0x00, 0xed, 0x32, 0xb2, 0x07,
    0x06, 0x01, 0xfc, 0x3c, 0xb2, 0x01, 0x02, 0x00, 0xed, 0x32, 0x33, 0x11, 0x21, 0x11, 0x25, 0x21, 0x11, 0x21, 0x3b,
    0x01, 0xda, 0xfe, 0x61, 0x01, 0x64, 0xfe, 0x9c, 0x04, 0xaa, 0xfb, 0x56, 0x3b, 0x04, 0x34, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0xff, 0x72, 0x05, 0x7f, 0x05, 0x8e, 0x00, 0x0b, 0x00, 0x00, 0x09, 0x01, 0x06, 0x26, 0x35, 0x11, 0x34,
    0x36, 0x17, 0x01, 0x16, 0x14, 0x05, 0x68, 0xfa, 0xd0, 0x17, 0x21, 0x21, 0x17, 0x05, 0x30, 0x17, 0x02, 0x61, 0xfd,
    0x1e, 0x0d, 0x14, 0x1a, 0x05, 0xc0, 0x1a, 0x14, 0x0d, 0xfd, 0x1e, 0x0d, 0x24, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00,
    0x00, 0xff, 0x80, 0x06, 0x00, 0x05, 0x80, 0x00, 0x0f, 0x00, 0x1f, 0x00, 0x00, 0x01, 0x11, 0x14, 0x06, 0x23, 0x21,
    0x22, 0x26, 0x35, 0x11, 0x34, 0x36, 0x33, 0x21, 0x32, 0x16, 0x05, 0x11, 0x14, 0x06, 0x23, 0x21, 0x22, 0x26, 0x35,
    0x11, 0x34, 0x36, 0x33, 0x21, 0x32, 0x16, 0x06, 0x00, 0x26, 0x1a, 0xfe, 0x00, 0x1a, 0x26, 0x26, 0x1a, 0x02, 0x00,
    0x1a, 0x26, 0xfc, 0x80, 0x26, 0x1a, 0xfe, 0x00, 0x1a, 0x26, 0x26, 0x1a, 0x02, 0x00, 0x1a, 0x26, 0x05, 0x40, 0xfa,
    0x80, 0x1a, 0x26, 0x26, 0x1a, 0x05, 0x80, 0x1a, 0x26, 0x26, 0x1a, 0xfa, 0x80, 0x1a, 0x26, 0x26, 0x1a, 0x05, 0x80,
    0x1a, 0x26, 0x26, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0xff, 0x80, 0x06, 0x00, 0x05, 0x80, 0x00, 0x0f, 0x00,
    0x00, 0x01, 0x11, 0x14, 0x06, 0x23, 0x21, 0x22, 0x26, 0x35, 0x11, 0x34, 0x36, 0x33, 0x21, 0x32, 0x16, 0x06, 0x00,
    0x26, 0x1a, 0xfa, 0x80, 0x1a, 0x26, 0x26, 0x1a, 0x05, 0x80, 0x1a, 0x26, 0x05, 0x40, 0xfa, 0x80, 0x1a, 0x26, 0x26,
    0x1a, 0x05, 0x80, 0x1a, 0x26, 0x26, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x05, 0x80, 0x05, 0x80,
    0x00, 0x23, 0x00, 0x00, 0x01, 0x15, 0x14, 0x06, 0x23, 0x21, 0x11, 0x14, 0x06, 0x2b, 0x01, 0x22, 0x26, 0x35, 0x11,
    0x21, 0x22, 0x26, 0x3d, 0x01, 0x34, 0x36, 0x33, 0x21, 0x11, 0x34, 0x36, 0x3b, 0x01, 0x32, 0x16, 0x15, 0x11, 0x21,
    0x32, 0x16, 0x05, 0x80, 0x38, 0x28, 0xfe, 0x60, 0x38, 0x28, 0xc0, 0x28, 0x38, 0xfe, 0x60, 0x28, 0x38, 0x38, 0x28,
    0x01, 0xa0, 0x38, 0x28, 0xc0, 0x28, 0x38, 0x01, 0xa0, 0x28, 0x38, 0x03, 0x20, 0xc0, 0x28, 0x38, 0xfe, 0x60, 0x28,
    0x38, 0x38, 0x28, 0x01, 0xa0, 0x38, 0x28, 0xc0, 0x28, 0x38, 0x01, 0xa0, 0x28, 0x38, 0x38, 0x28, 0xfe, 0x60, 0x38,
    0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x02, 0x00, 0x05, 0x80, 0x03, 0x80, 0x00, 0x0f, 0x00, 0x00, 0x01, 0x15,
    0x14, 0x06, 0x23, 0x21, 0x22, 0x26, 0x3d, 0x01, 0x34, 0x36, 0x33, 0x21, 0x32, 0x16, 0x05, 0x80, 0x38, 0x28, 0xfb,
    0x40, 0x28, 0x38, 0x38, 0x28, 0x04, 0xc0, 0x28, 0x38, 0x03, 0x20, 0xc0, 0x28, 0x38, 0x38, 0x28, 0xc0, 0x28, 0x38,
    0x38, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x07, 0x00, 0x04, 0x80, 0x00, 0x11, 0x00, 0x21, 0x00, 0x31, 0x00,
    0x00, 0x01, 0x26, 0x27, 0x16, 0x15, 0x14, 0x00, 0x20, 0x00, 0x35, 0x34, 0x37, 0x06, 0x07, 0x16, 0x04, 0x20, 0x24,
    0x00, 0x34, 0x26, 0x23, 0x22, 0x06, 0x15, 0x14, 0x16, 0x32, 0x36, 0x35, 0x34, 0x36, 0x33, 0x32, 0x00, 0x14, 0x07,
    0x06, 0x00, 0x20, 0x00, 0x27, 0x26, 0x34, 0x37, 0x36, 0x00, 0x20, 0x00, 0x17, 0x06, 0x80, 0x98, 0xe5, 0x3d, 0xfe,
    0xf9, 0xfe, 0x8e, 0xfe, 0xf9, 0x3d, 0xe5, 0x98, 0x85, 0x01, 0x91, 0x01, 0xd4, 0x01, 0x91, 0xfd, 0xb5, 0x1c, 0x14,
    0x7d, 0xb3, 0x1c, 0x28, 0x1c, 0x7a, 0x56, 0x14, 0x03, 0x6c, 0x14, 0x8c, 0xfe, 0x27, 0xfd, 0xf2, 0xfe, 0x27, 0x8c,
    0x14, 0x14, 0x8c, 0x01, 0xd9, 0x02, 0x0e, 0x01, 0xd9, 0x8c, 0x02, 0x40, 0xec, 0x75, 0x68, 0x79, 0xb9, 0xfe, 0xf9,
    0x01, 0x07, 0xb9, 0x79, 0x68, 0x75, 0xec, 0xcd, 0xf3, 0xf3, 0x02, 0x39, 0x28, 0x1c, 0xb3, 0x7d, 0x14, 0x1c, 0x1c,
    0x14, 0x56, 0x7a, 0xfe, 0xd2, 0x44, 0x23, 0xe6, 0xfe, 0xeb, 0x01, 0x16, 0xe5, 0x23, 0x44, 0x23, 0xe5, 0x01, 0x16,
    0xfe, 0xea, 0xe5, 0x00, 0x03, 0x00, 0x00, 0xff, 0x8f, 0x07, 0x26, 0x04, 0xf1, 0x00, 0x14, 0x00, 0x24, 0x00, 0x39,
    0x00, 0x00, 0x25, 0x07, 0x06, 0x22, 0x27, 0x01, 0x26, 0x34, 0x37, 0x01, 0x36, 0x32, 0x1f, 0x01, 0x16, 0x14, 0x07,
    0x09, 0x01, 0x16, 0x14, 0x09, 0x01, 0x0e, 0x01, 0x2f, 0x01, 0x2e, 0x01, 0x37, 0x01, 0x3e, 0x01, 0x1f, 0x01, 0x1e,
    0x01, 0x09, 0x01, 0x06, 0x22, 0x2f, 0x01, 0x26, 0x34, 0x37, 0x09, 0x01, 0x26, 0x34, 0x3f, 0x01, 0x36, 0x32, 0x17,
    0x01, 0x16, 0x14, 0x02, 0x3c, 0x32, 0x0a, 0x1a, 0x0a, 0xfe, 0x2e, 0x0a, 0x0a, 0x01, 0xd2, 0x0a, 0x1a, 0x0a, 0x32,
    0x0a, 0x0a, 0xfe, 0x77, 0x01, 0x89, 0x0a, 0x02, 0x45, 0xfe, 0x8b, 0x04, 0x17, 0x0c, 0x3e, 0x0d, 0x0d, 0x04, 0x01,
    0x75, 0x04, 0x17, 0x0c, 0x3e, 0x0d, 0x0d, 0x02, 0x8d, 0xfe, 0x2e, 0x0a, 0x1a, 0x0a, 0x32, 0x0a, 0x0a, 0x01, 0x89,
    0xfe, 0x77, 0x0a, 0x0a, 0x32, 0x0a, 0x1a, 0x0a, 0x01, 0xd2, 0x0a, 0x89, 0x32, 0x0a, 0x0a, 0x01, 0xd2, 0x0a, 0x1a,
    0x0a, 0x01, 0xd2, 0x0a, 0x0a, 0x32, 0x0a, 0x1a, 0x0a, 0xfe, 0x77, 0xfe, 0x77, 0x0a, 0x1a, 0x04, 0x21, 0xfa, 0xf5,
    0x0d, 0x0d, 0x04, 0x11, 0x04, 0x17, 0x0d, 0x05, 0x0b, 0x0d, 0x0d, 0x04, 0x11, 0x04, 0x17, 0xfd, 0x68, 0xfe, 0x2e,
    0x0a, 0x0a, 0x32, 0x0a, 0x1a, 0x0a, 0x01, 0x89, 0x01, 0x89, 0x0a, 0x1a, 0x0a, 0x32, 0x0a, 0x0a, 0xfe, 0x2e, 0x0a,
    0x1a, 0x00, 0x00, 0x00, 0x00, 0x0e, 0x00, 0xae, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xd2, 0x01,
    0xa6, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x0b, 0x02, 0x91, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x02, 0x00, 0x07, 0x02, 0xad, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x27, 0x03, 0x05, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x0b, 0x03, 0x45, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05,
    0x00, 0x10, 0x03, 0x73, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x00, 0x0b, 0x03, 0x9c, 0x00, 0x03, 0x00,
    0x01, 0x04, 0x09, 0x00, 0x00, 0x01, 0xa4, 0x00, 0x00, 0x00, 0x03, 0x00, 0x01, 0x04, 0x09, 0x00, 0x01, 0x00, 0x16,
    0x02, 0x79, 0x00, 0x03, 0x00, 0x01, 0x04, 0x09, 0x00, 0x02, 0x00, 0x0e, 0x02, 0x9d, 0x00, 0x03, 0x00, 0x01, 0x04,
    0x09, 0x00, 0x03, 0x00, 0x4e, 0x02, 0xb5, 0x00, 0x03, 0x00, 0x01, 0x04, 0x09, 0x00, 0x04, 0x00, 0x16, 0x03, 0x2d,
    0x00, 0x03, 0x00, 0x01, 0x04, 0x09, 0x00, 0x05, 0x00, 0x20, 0x03, 0x51, 0x00, 0x03, 0x00, 0x01, 0x04, 0x09, 0x00,
    0x06, 0x00, 0x16, 0x03, 0x84, 0x00, 0x54, 0x00, 0x68, 0x00, 0x65, 0x00, 0x20, 0x00, 0x46, 0x00, 0x6f, 0x00, 0x72,
    0x00, 0x6b, 0x00, 0x20, 0x00, 0x41, 0x00, 0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00,
    0x20, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x6e, 0x00, 0x74, 0x00, 0x20, 0x00, 0x69, 0x00, 0x73, 0x00, 0x20, 0x00, 0x6c,
    0x00, 0x69, 0x00, 0x63, 0x00, 0x65, 0x00, 0x6e, 0x00, 0x73, 0x00, 0x65, 0x00, 0x64, 0x00, 0x20, 0x00, 0x75, 0x00,
    0x6e, 0x00, 0x64, 0x00, 0x65, 0x00, 0x72, 0x00, 0x20, 0x00, 0x74, 0x00, 0x68, 0x00, 0x65, 0x00, 0x20, 0x00, 0x53,
    0x00, 0x49, 0x00, 0x4c, 0x00, 0x20, 0x00, 0x4f, 0x00, 0x46, 0x00, 0x4c, 0x00, 0x20, 0x00, 0x31, 0x00, 0x2e, 0x00,
    0x31, 0x00, 0x20, 0x00, 0x28, 0x00, 0x68, 0x00, 0x74, 0x00, 0x74, 0x00, 0x70, 0x00, 0x3a, 0x00, 0x2f, 0x00, 0x2f,
    0x00, 0x73, 0x00, 0x63, 0x00, 0x72, 0x00, 0x69, 0x00, 0x70, 0x00, 0x74, 0x00, 0x73, 0x00, 0x2e, 0x00, 0x73, 0x00,
    0x69, 0x00, 0x6c, 0x00, 0x2e, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x67, 0x00, 0x2f, 0x00, 0x4f, 0x00, 0x46, 0x00, 0x4c,
    0x00, 0x29, 0x00, 0x2e, 0x00, 0x20, 0x00, 0x46, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b, 0x00, 0x20, 0x00, 0x41, 0x00,
    0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x20, 0x00, 0x69, 0x00, 0x73, 0x00, 0x20,
    0x00, 0x61, 0x00, 0x20, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b, 0x00, 0x20, 0x00, 0x62, 0x00, 0x61, 0x00,
    0x73, 0x00, 0x65, 0x00, 0x64, 0x00, 0x20, 0x00, 0x6f, 0x00, 0x66, 0x00, 0x20, 0x00, 0x6f, 0x00, 0x66, 0x00, 0x66,
    0x00, 0x20, 0x00, 0x46, 0x00, 0x6f, 0x00, 0x6e, 0x00, 0x74, 0x00, 0x20, 0x00, 0x41, 0x00, 0x77, 0x00, 0x65, 0x00,
    0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x20, 0x00, 0x34, 0x00, 0x2e, 0x00, 0x37, 0x00, 0x2e, 0x00, 0x30,
    0x00, 0x20, 0x00, 0x62, 0x00, 0x79, 0x00, 0x20, 0x00, 0x44, 0x00, 0x61, 0x00, 0x76, 0x00, 0x65, 0x00, 0x20, 0x00,
    0x47, 0x00, 0x61, 0x00, 0x6e, 0x00, 0x64, 0x00, 0x79, 0x00, 0x2e, 0x00, 0x20, 0x00, 0x4d, 0x00, 0x6f, 0x00, 0x72,
    0x00, 0x65, 0x00, 0x20, 0x00, 0x69, 0x00, 0x6e, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x20, 0x00, 0x6f, 0x00, 0x6e, 0x00,
    0x20, 0x00, 0x6c, 0x00, 0x69, 0x00, 0x63, 0x00, 0x65, 0x00, 0x6e, 0x00, 0x73, 0x00, 0x65, 0x00, 0x73, 0x00, 0x20,
    0x00, 0x61, 0x00, 0x74, 0x00, 0x20, 0x00, 0x68, 0x00, 0x74, 0x00, 0x74, 0x00, 0x70, 0x00, 0x73, 0x00, 0x3a, 0x00,
    0x2f, 0x00, 0x2f, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b, 0x00, 0x61, 0x00, 0x77, 0x00, 0x65, 0x00, 0x73,
    0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x2e, 0x00, 0x67, 0x00, 0x69, 0x00, 0x74, 0x00, 0x68, 0x00, 0x75, 0x00,
    0x62, 0x00, 0x2e, 0x00, 0x69, 0x00, 0x6f, 0x00, 0x00, 0x54, 0x68, 0x65, 0x20, 0x46, 0x6f, 0x72, 0x6b, 0x20, 0x41,
    0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x20, 0x66, 0x6f, 0x6e, 0x74, 0x20, 0x69, 0x73, 0x20, 0x6c, 0x69, 0x63, 0x65,
    0x6e, 0x73, 0x65, 0x64, 0x20, 0x75, 0x6e, 0x64, 0x65, 0x72, 0x20, 0x74, 0x68, 0x65, 0x20, 0x53, 0x49, 0x4c, 0x20,
    0x4f, 0x46, 0x4c, 0x20, 0x31, 0x2e, 0x31, 0x20, 0x28, 0x68, 0x74, 0x74, 0x70, 0x3a, 0x2f, 0x2f, 0x73, 0x63, 0x72,
    0x69, 0x70, 0x74, 0x73, 0x2e, 0x73, 0x69, 0x6c, 0x2e, 0x6f, 0x72, 0x67, 0x2f, 0x4f, 0x46, 0x4c, 0x29, 0x2e, 0x20,
    0x46, 0x6f, 0x72, 0x6b, 0x20, 0x41, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x20, 0x69, 0x73, 0x20, 0x61, 0x20, 0x66,
    0x6f, 0x72, 0x6b, 0x20, 0x62, 0x61, 0x73, 0x65, 0x64, 0x20, 0x6f, 0x66, 0x20, 0x6f, 0x66, 0x66, 0x20, 0x46, 0x6f,
    0x6e, 0x74, 0x20, 0x41, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x20, 0x34, 0x2e, 0x37, 0x2e, 0x30, 0x20, 0x62, 0x79,
    0x20, 0x44, 0x61, 0x76, 0x65, 0x20, 0x47, 0x61, 0x6e, 0x64, 0x79, 0x2e, 0x20, 0x4d, 0x6f, 0x72, 0x65, 0x20, 0x69,
    0x6e, 0x66, 0x6f, 0x20, 0x6f, 0x6e, 0x20, 0x6c, 0x69, 0x63, 0x65, 0x6e, 0x73, 0x65, 0x73, 0x20, 0x61, 0x74, 0x20,
    0x68, 0x74, 0x74, 0x70, 0x73, 0x3a, 0x2f, 0x2f, 0x66, 0x6f, 0x72, 0x6b, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65,
    0x2e, 0x67, 0x69, 0x74, 0x68, 0x75, 0x62, 0x2e, 0x69, 0x6f, 0x00, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b,
    0x00, 0x61, 0x00, 0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x00, 0x66, 0x6f, 0x72,
    0x6b, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x00, 0x00, 0x52, 0x00, 0x65, 0x00, 0x67, 0x00, 0x75, 0x00, 0x6c,
    0x00, 0x61, 0x00, 0x72, 0x00, 0x00, 0x52, 0x65, 0x67, 0x75, 0x6c, 0x61, 0x72, 0x00, 0x00, 0x46, 0x00, 0x6f, 0x00,
    0x6e, 0x00, 0x74, 0x00, 0x46, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x67, 0x00, 0x65, 0x00, 0x20, 0x00, 0x32, 0x00, 0x2e,
    0x00, 0x30, 0x00, 0x20, 0x00, 0x3a, 0x00, 0x20, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b, 0x00, 0x61, 0x00,
    0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x20, 0x00, 0x3a, 0x00, 0x20, 0x00, 0x32,
    0x00, 0x34, 0x00, 0x2d, 0x00, 0x36, 0x00, 0x2d, 0x00, 0x32, 0x00, 0x30, 0x00, 0x32, 0x00, 0x30, 0x00, 0x00, 0x46,
    0x6f, 0x6e, 0x74, 0x46, 0x6f, 0x72, 0x67, 0x65, 0x20, 0x32, 0x2e, 0x30, 0x20, 0x3a, 0x20, 0x66, 0x6f, 0x72, 0x6b,
    0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x20, 0x3a, 0x20, 0x32, 0x34, 0x2d, 0x36, 0x2d, 0x32, 0x30, 0x32, 0x30,
    0x00, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00, 0x6b, 0x00, 0x61, 0x00, 0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f,
    0x00, 0x6d, 0x00, 0x65, 0x00, 0x00, 0x66, 0x6f, 0x72, 0x6b, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x00, 0x00,
    0x56, 0x00, 0x65, 0x00, 0x72, 0x00, 0x73, 0x00, 0x69, 0x00, 0x6f, 0x00, 0x6e, 0x00, 0x20, 0x00, 0x30, 0x00, 0x30,
    0x00, 0x31, 0x00, 0x2e, 0x00, 0x30, 0x00, 0x30, 0x00, 0x30, 0x00, 0x20, 0x00, 0x00, 0x56, 0x65, 0x72, 0x73, 0x69,
    0x6f, 0x6e, 0x20, 0x30, 0x30, 0x31, 0x2e, 0x30, 0x30, 0x30, 0x20, 0x00, 0x00, 0x66, 0x00, 0x6f, 0x00, 0x72, 0x00,
    0x6b, 0x00, 0x61, 0x00, 0x77, 0x00, 0x65, 0x00, 0x73, 0x00, 0x6f, 0x00, 0x6d, 0x00, 0x65, 0x00, 0x00, 0x66, 0x6f,
    0x72, 0x6b, 0x61, 0x77, 0x65, 0x73, 0x6f, 0x6d, 0x65, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x7c, 0x00, 0x59, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x0b, 0x00, 0x00, 0x00, 0x01, 0x00, 0x02, 0x00, 0x03, 0x01, 0x02, 0x01,
    0x03, 0x01, 0x04, 0x01, 0x05, 0x01, 0x06, 0x01, 0x07, 0x01, 0x08, 0x07, 0x75, 0x6e, 0x69, 0x46, 0x30, 0x30, 0x30,
    0x07, 0x75, 0x6e, 0x69, 0x46, 0x30, 0x30, 0x31, 0x07, 0x75, 0x6e, 0x69, 0x46, 0x30, 0x30, 0x32, 0x07, 0x75, 0x6e,
    0x69, 0x46, 0x30, 0x30, 0x33, 0x07, 0x75, 0x6e, 0x69, 0x46, 0x30, 0x30, 0x34, 0x07, 0x75, 0x6e, 0x69, 0x46, 0x30,
    0x30, 0x35, 0x07, 0x75, 0x6e, 0x69, 0x46, 0x30, 0x30, 0x36, 0x00, 0x00, 0x00, 0x01, 0xff, 0xff, 0x00, 0x02, 0x00,
    0x01, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x16, 0x00, 0x00, 0x00, 0x02, 0x00, 0x01, 0x00, 0x01, 0x00, 0x0a,
    0x00, 0x01, 0x00, 0x04, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0x00, 0xdb, 0xcc, 0xbf, 0x7d, 0x00, 0x00, 0x00, 0x00, 0xdb, 0x18, 0xe0, 0x91, 0x00, 0x00, 0x00, 0x00, 0xdc, 0x9b,
    0x30, 0x5d,
};
