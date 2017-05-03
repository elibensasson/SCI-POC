#include <gtest/gtest.h>
#include <gadgetlib/gadget.hpp>
#include <gadgetlib/protoboard.hpp>

#define ROUNDS 111
#define REGISTER_SIZE 16

using namespace gadgetlib;


namespace{
	void setValueToUnpacked(ProtoboardPtr pb, Algebra::UnpackedWord unpacked, int value){
		for (int i = 0; i < unpacked.size(); i++){
			pb->val(unpacked[i]) = value & 1u ? Algebra::one() : Algebra::zero();
			value >>= 1;
		}
	}
	/*
	TEST(gadget, greaterEqual_gadget){
		ProtoboardPtr pb = Protoboard::create();
		Algebra::UnpackedWord a(REGISTER_SIZE, "a");
		Algebra::UnpackedWord b(REGISTER_SIZE, "b");
		Algebra::FlagVariable flagGreater;
		Algebra::FlagVariable flagEqual;
		GadgetPtr greaterGadget = GreaterEqual_Gadget::create(pb,a,b,flagGreater,flagEqual,Opcode::NONE);
		greaterGadget->generateConstraints();
		for (int i = 0; i < ROUNDS; i++){
			for (int j = 0; j < ROUNDS; j++){
				setValueToUnpacked(pb, a, i);
				setValueToUnpacked(pb, b, j);
				greaterGadget->generateWitness();
				if (i > j){
					EXPECT_EQ(pb->val(flagGreater), Algebra::one());
					EXPECT_EQ(pb->val(flagEqual), Algebra::zero());
				}
				else{
					if (i < j){
						EXPECT_EQ(pb->val(flagGreater), Algebra::zero());
						EXPECT_EQ(pb->val(flagEqual), Algebra::zero());
					}
					else{ // i == j
						EXPECT_EQ(pb->val(flagGreater), Algebra::zero());
						EXPECT_EQ(pb->val(flagEqual), Algebra::one());
					}
				}
				EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
			}
		}
	}
	*/
	TEST(gadget, multiplicationPacking){
		ProtoboardPtr pb = Protoboard::create();
		Algebra::UnpackedWord a_unpacked(REGISTER_SIZE, "a_unpacked");
		Algebra::UnpackedWord a_results(REGISTER_SIZE, "a_results");

		GadgetPtr packingGadget = MultiplicationPacking_Gadget::create(pb, a_unpacked, a_results,
			a_results[a_results.size()-1], false, PackingMode::PACK, Opcode::NONE);
		packingGadget->generateConstraints();
		const Algebra::FElem g = Algebra::FElem(getGF2E_X());
		Algebra::FElem x_i = Algebra::one();
		for (int i = 0; i < ROUNDS; i++){
			setValueToUnpacked(pb, a_unpacked, i);
			std::dynamic_pointer_cast<MultiplicationPacking_Gadget>(packingGadget)->generateWitness();
			EXPECT_EQ(pb->val(a_results[a_results.size()-1]), x_i);
			EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
			x_i *= g;
		}
		GadgetPtr unpackingGadget = MultiplicationPacking_Gadget::create(pb, a_unpacked, a_results,
			a_results[a_results.size()-1], false, PackingMode::UNPACK, Opcode::NONE);
		unpackingGadget->generateConstraints();
		x_i = Algebra::one();
		for (int i = 0; i < ROUNDS; i++){
			pb->val(a_results[a_results.size()-1]) = x_i;
			std::dynamic_pointer_cast<MultiplicationPacking_Gadget>(unpackingGadget)->generateWitness(i);
			int value = i;
			for (int j = 0; j < a_unpacked.size(); j++){
				Algebra::FElem val = (value & 1u ? Algebra::one() : Algebra::zero());
				EXPECT_EQ(pb->val(a_unpacked[j]), val);
				value >>= 1u;
			}
			EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
			x_i *= g;
		}
	}
	TEST(gadget, mult){
		ProtoboardPtr pb = Protoboard::create();
		Algebra::UnpackedWord input1(REGISTER_SIZE, "inp1");
		Algebra::UnpackedWord input2(REGISTER_SIZE, "inp2");
		Algebra::UnpackedWord partials1(REGISTER_SIZE, "prt1");
		Algebra::UnpackedWord partials2(REGISTER_SIZE, "prt2");
		Algebra::UnpackedWord result(REGISTER_SIZE, "res");
		Algebra::UnpackedWord packed_results(REGISTER_SIZE, "packed_res");

		GadgetPtr packGadget = MultiplicationPacking_Gadget::create(pb, result, packed_results,
			packed_results[packed_results.size()-1], false, PackingMode::PACK, Opcode::NONE);
		GadgetPtr multGadget = Multiplication_Gadget::create(pb, input1, input2, partials1, partials2, false, Opcode::NONE);
		multGadget->generateConstraints();
		packGadget->generateConstraints();
		for (int i = 0; i < ROUNDS; i++){
			setValueToUnpacked(pb, input1, i);
			for (int j = 0; j < ROUNDS; j++){
				setValueToUnpacked(pb, input2, j);
				setValueToUnpacked(pb, result, i*j);
				multGadget->generateWitness();
				std::dynamic_pointer_cast<MultiplicationPacking_Gadget>(packGadget)->generateWitness();
				EXPECT_EQ(pb->val(packed_results[packed_results.size()-1]), pb->val(partials2[0]));
				EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
			}
		}
	}
}