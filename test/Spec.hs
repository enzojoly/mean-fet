-- |
-- Module      : Spec
-- Description : Test suite entry point.

import Test.Tasty (defaultMain, testGroup)

import qualified UnitSpec
import qualified PropertySpec
import qualified RegressionSpec
import qualified ConsistencySpec
import qualified BoundarySpec
import qualified EncounterSpec
import qualified GridSpec
import qualified ReferenceSpec
import qualified EstimateSpec
import qualified OracleSpec
import qualified DefectSpec
import qualified AbsorbingSpec
import qualified CellsSpec
import qualified PrimitiveSpec
import qualified DecompositionSpec
import qualified SensitivitySpec
import qualified SlopeSpec

main :: IO ()
main = defaultMain $ testGroup "Bivia"
    [ UnitSpec.tests
    , PropertySpec.tests
    , RegressionSpec.tests
    , ConsistencySpec.tests
    , BoundarySpec.tests
    , EncounterSpec.tests
    , GridSpec.tests
    , ReferenceSpec.tests
    , EstimateSpec.tests
    , OracleSpec.tests
    , DefectSpec.tests
    , AbsorbingSpec.tests
    , CellsSpec.tests
    , PrimitiveSpec.tests
    , DecompositionSpec.tests
    , SensitivitySpec.tests
    , SlopeSpec.tests
    ]
