################################
# Complex Operator Definitions #
################################
@present ExtendedOperators <: ExtendedFreeExtCalc2D begin
  X::Space
  neg::Hom(DualForm1(X), DualForm1(X)) # negative
  neg₁::Hom(Form1(X), Form1(X)) # negates 1-forms
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # negates 1-forms
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  i₀::Hom(Form1(X)⊗DualForm2(X), DualForm1(X))
  i₁::Hom(Form1(X)⊗DualForm1(X), DualForm0(X))
    
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  i₀′::Hom(Form1(X)⊗Form2(X), Form1(X))
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  ∧₁₁′::Hom(Form1(X)⊗DualForm1(X), DualForm2(X))
  ∧₁₀′::Hom(Form1(X)⊗DualForm0(X), Form1(X))
end

rules = gen_dec_rules()

# Lie derivative between two primal 1-forms
Lie1Imp′ = @decapode ExtendedOperators begin
  (F1, F1′, Ḟ1)::Form1{X}
  Ḟ1 == i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
end
lie1_imp′ = diag2dwd(Lie1Imp′, in_vars = [:F1, :F1′], out_vars = [:Ḟ1])
rules[:L₁′] = lie1_imp′

# Internal product between a primal 1-form and a primal 2-form
I0Imp′ = @decapode ExtendedOperators begin
  (F1, Ḟ1)::Form1{X}
  F2::Form2{X}
  Ḟ1 == neg₁(∧₁₀′(F1, ⋆₂{X}(F2)))
end
i0_imp′ = diag2dwd(I0Imp′, in_vars = [:F1, :F2], out_vars = [:Ḟ1])
rules[:i₀′] = i0_imp′

# Internal product between two primal 1-forms
I1Imp′ = @decapode ExtendedOperators begin
  (F1, F1′)::Form1{X}
  F0::Form0{X}
  F0 == ⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
end
i1_imp′ = diag2dwd(I1Imp′, in_vars = [:F1, :F1′], out_vars = [:F0])
rules[:i₁′] = i1_imp′


