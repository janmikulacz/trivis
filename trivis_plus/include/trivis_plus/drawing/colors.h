/**
 * File:    colors.h
 *
 * Date:    25.08.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DRAWING_COLORS_H_
#define TRIVIS_PLUS_DRAWING_COLORS_H_

namespace trivis_plus::drawing {

struct RGB {
  int r = 255;
  int g = 255;
  int b = 255;
};

/// ##########
/// ## Red ##
/// #########

static constexpr RGB kColorLightSalmon = {255, 160, 122};
static constexpr RGB kColorSalmon = {250, 128, 114};
static constexpr RGB kColorDarkSalmon = {233, 150, 122};
static constexpr RGB kColorLightCoral = {240, 128, 128};
static constexpr RGB kColorIndianRed = {205, 92, 92};
static constexpr RGB kColorCrimson = {220, 20, 60};
static constexpr RGB kColorFirebrick = {178, 34, 34};
static constexpr RGB kColorRed = {255, 0, 0};
static constexpr RGB kColorDarkRed = {139, 0, 0};

/// ############
/// ## Orange ##
/// ############

static constexpr RGB kColorCoral = {255, 127, 80};
static constexpr RGB kColorTomato = {255, 99, 71};
static constexpr RGB kColorOrangeRed = {255, 69, 0};
static constexpr RGB kColorGold = {255, 215, 0};
static constexpr RGB kColorOrange = {255, 165, 0};
static constexpr RGB kColorDarkOrange = {255, 140, 0};

/// ############
/// ## Yellow ##
/// ############

static constexpr RGB kColorLightYellow = {255, 255, 224};
static constexpr RGB kColorLemonChiffon = {255, 250, 205};
static constexpr RGB kColorLightGoldenRodYellow = {250, 250, 210};
static constexpr RGB kColorPapayaWhip = {255, 239, 213};
static constexpr RGB kColorMoccasin = {255, 228, 181};
static constexpr RGB kColorPeachPuff = {255, 218, 185};
static constexpr RGB kColorPaleGoldenRod = {238, 232, 170};
static constexpr RGB kColorKhaki = {240, 230, 140};
static constexpr RGB kColorDarkKhaki = {189, 183, 107};
static constexpr RGB kColorYellow = {255, 255, 0};

/// ###########
/// ## Green ##
/// ###########

static constexpr RGB kColorLawnGreen = {124, 252, 0};
static constexpr RGB kColorChartreuse = {127, 255, 0};
static constexpr RGB kColorLimeGreen = {50, 205, 50};
static constexpr RGB kColorLime = {0, 255, 0};
static constexpr RGB kColorForestGreen = {34, 139, 34};
static constexpr RGB kColorGreen = {0, 128, 0};
static constexpr RGB kColorDarkGreen = {0, 100, 0};
static constexpr RGB kColorGreenYellow = {173, 255, 47};
static constexpr RGB kColorYellowGreen = {154, 205, 50};
static constexpr RGB kColorSpringGreen = {0, 255, 127};
static constexpr RGB kColorMediumSpringGreen = {0, 250, 154};
static constexpr RGB kColorLightGreen = {144, 238, 144};
static constexpr RGB kColorPaleGreen = {152, 251, 152};
static constexpr RGB kColorDarkSeaGreen = {143, 188, 143};
static constexpr RGB kColorMediumSeaGreen = {60, 179, 113};
static constexpr RGB kColorSeaGreen = {46, 139, 87};
static constexpr RGB kColorOlive = {128, 128, 0};
static constexpr RGB kColorDarkOliveGreen = {85, 107, 47};
static constexpr RGB kColorOliveDrab = {107, 142, 35};

/// ##########
/// ## Cyan ##
/// ##########

static constexpr RGB kColorLightCyan = {224, 255, 255};
static constexpr RGB kColorCyan = {0, 255, 255};
static constexpr RGB kColorAqua = {0, 255, 255};
static constexpr RGB kColorAquamarine = {127, 255, 212};
static constexpr RGB kColorMediumAquamarine = {102, 205, 170};
static constexpr RGB kColorPaleTurquoise = {175, 238, 238};
static constexpr RGB kColorTurquoise = {64, 224, 208};
static constexpr RGB kColorMediumTurquoise = {72, 209, 204};
static constexpr RGB kColorDarkTurquoise = {0, 206, 209};
static constexpr RGB kColorLightSeaGreen = {32, 178, 170};
static constexpr RGB kColorCadetBlue = {95, 158, 160};
static constexpr RGB kColorDarkCyan = {0, 139, 139};
static constexpr RGB kColorTeal = {0, 128, 128};

/// ##########
/// ## Blue ##
/// ##########

static constexpr RGB kColorPowderBlue = {176, 224, 230};
static constexpr RGB kColorLightBlue = {173, 216, 230};
static constexpr RGB kColorLightSkyBlue = {135, 206, 250};
static constexpr RGB kColorSkyBlue = {135, 206, 235};
static constexpr RGB kColorDeepSkyBlue = {0, 191, 255};
static constexpr RGB kColorLightSteelBlue = {176, 196, 222};
static constexpr RGB kColorDodgerBlue = {30, 144, 255};
static constexpr RGB kColorCornFlowerBlue = {100, 149, 237};
static constexpr RGB kColorSteelBlue = {70, 130, 180};
static constexpr RGB kColorRoyalBlue = {65, 105, 225};
static constexpr RGB kColorBlue = {0, 0, 255};
static constexpr RGB kColorMediumBlue = {0, 0, 205};
static constexpr RGB kColorDarkBlue = {0, 0, 139};
static constexpr RGB kColorNavy = {0, 0, 128};
static constexpr RGB kColorMidnightBlue = {25, 25, 112};
static constexpr RGB kColorMediumSlateBlue = {123, 104, 238};
static constexpr RGB kColorSlateBlue = {106, 90, 205};
static constexpr RGB kColorDarkSlateBlue = {72, 61, 139};

/// ############
/// ## Purple ##
/// ############

static constexpr RGB kColorLavender = {230, 230, 250};
static constexpr RGB kColorThistle = {216, 191, 216};
static constexpr RGB kColorPlum = {221, 160, 221};
static constexpr RGB kColorViolet = {238, 130, 238};
static constexpr RGB kColorOrchid = {218, 112, 214};
static constexpr RGB kColorFuchsia = {255, 0, 255};
static constexpr RGB kColorMagenta = {255, 0, 255};
static constexpr RGB kColorMediumOrchid = {186, 85, 211};
static constexpr RGB kColorMediumPurple = {147, 112, 219};
static constexpr RGB kColorBlueViolet = {138, 43, 226};
static constexpr RGB kColorDarkViolet = {148, 0, 211};
static constexpr RGB kColorDarkOrchid = {153, 50, 204};
static constexpr RGB kColorDarkMagenta = {139, 0, 139};
static constexpr RGB kColorPurple = {128, 0, 128};
static constexpr RGB kColorIndigo = {75, 0, 130};

/// ##########
/// ## Pink ##
/// ##########

static constexpr RGB kColorPink = {255, 192, 203};
static constexpr RGB kColorLightPink = {255, 182, 193};
static constexpr RGB kColorHotPink = {255, 105, 180};
static constexpr RGB kColorDeepPink = {255, 20, 147};
static constexpr RGB kColorPaleVioletRed = {219, 112, 147};
static constexpr RGB kColorMediumVioletRed = {199, 21, 133};

/// ###########
/// ## White ##
/// ###########

static constexpr RGB kColorWhite = {255, 255, 255};
static constexpr RGB kColorSnow = {255, 250, 250};
static constexpr RGB kColorHoneyDew = {240, 255, 240};
static constexpr RGB kColorMintCream = {245, 255, 250};
static constexpr RGB kColorAzure = {240, 255, 255};
static constexpr RGB kColorAliceBlue = {240, 248, 255};
static constexpr RGB kColorGhostWhite = {248, 248, 255};
static constexpr RGB kColorWhiteSmoke = {245, 245, 245};
static constexpr RGB kColorSeashell = {255, 245, 238};
static constexpr RGB kColorBeige = {245, 245, 220};
static constexpr RGB kColorOldLace = {253, 245, 230};
static constexpr RGB kColorFloralWhite = {255, 250, 240};
static constexpr RGB kColorIvory = {255, 255, 240};
static constexpr RGB kColorAntiqueWhite = {250, 235, 215};
static constexpr RGB kColorLinen = {250, 240, 230};
static constexpr RGB kColorLavenderBlush = {255, 240, 245};
static constexpr RGB kColorMistyRose = {255, 228, 225};

/// ##########
/// ## Gray ##
/// ##########

static constexpr RGB kColorGainsboro = {220, 220, 220};
static constexpr RGB kColorLightGray = {211, 211, 211};
static constexpr RGB kColorSilver = {192, 192, 192};
static constexpr RGB kColorDarkGray = {169, 169, 169};
static constexpr RGB kColorGray = {128, 128, 128};
static constexpr RGB kColorDimGray = {105, 105, 105};
static constexpr RGB kColorLightSlateGray = {119, 136, 153};
static constexpr RGB kColorSlateGray = {112, 128, 144};
static constexpr RGB kColorDarkSlateGray = {47, 79, 79};
static constexpr RGB kColorBlack = {0, 0, 0};

/// ###########
/// ## Brown ##
/// ###########

static constexpr RGB kColorCornsilk = {255, 248, 220};
static constexpr RGB kColorBlanchedAlmond = {255, 235, 205};
static constexpr RGB kColorBisque = {255, 228, 196};
static constexpr RGB kColorNavajoWhite = {255, 222, 173};
static constexpr RGB kColorWheat = {245, 222, 179};
static constexpr RGB kColorBurlywood = {222, 184, 135};
static constexpr RGB kColorTan = {210, 180, 140};
static constexpr RGB kColorRosyBrown = {188, 143, 143};
static constexpr RGB kColorSandyBrown = {244, 164, 96};
static constexpr RGB kColorGoldenrod = {218, 165, 32};
static constexpr RGB kColorPeru = {205, 133, 63};
static constexpr RGB kColorChocolate = {210, 105, 30};
static constexpr RGB kColorSaddleBrown = {139, 69, 19};
static constexpr RGB kColorSienna = {160, 82, 45};
static constexpr RGB kColorBrown = {165, 42, 42};
static constexpr RGB kColorMaroon = {128, 0, 0};

}

#endif //TRIVIS_PLUS_DRAWING_COLORS_H_
