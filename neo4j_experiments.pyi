from typing import Dict, Callable, TypeVar
import Bio
import combat_tb_neomodel.model.core

Feature = combat_tb_neomodel.model.core.Feature
SeqFeature = Bio.SeqFeature.SeqFeature
OptionalTransform = TypeVar('OptionalTransform', Callable([str], str), None)

def has_qualifier(feature: SeqFeature,
                  key: str,
                  value: str) -> bool: pass

def set_if_has(thing: Feature,
               feature: SeqFeature,
               key: str,
               transform: OptionalTransform) -> None: pass

def strip_id_colon(id: str) -> str: pass

def get_parent(feature: SeqFeature,
               thing: Feature,
               lookup_dict: Dict[str, Feature]) -> Feature: pass